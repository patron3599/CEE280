classdef CPBD_Analysis < handle

% This class contains the necesary methods to complete all computations
% involved in assesing the respone of the structure (i.e displacements, 
% suppport reactions, and element forces) to the loading it is subjected
% to. An instance of this class created when the funtion is called and
% contains the global structural properties that are nodes, elements,
% deflections, reactions, element forces, and flags. All these parameters
% will then be passed into the MASTAN program and subsequently retrieved to
% display the strucutre's response. 
    
    % Private properties go here
    properties (Access = private)

        eleobj
        element_nodes
        DEFL
        REACT
        ELE_FOR
        AFLAG
     
        %% Property Definitions 
        %  Note: Variables defined in ud_3d1el.m are not repeated here
        %
        %  eleobj  == The object of each element containing its respective 
        %             degrees of freedom and nodes. This object is then 
        %             passed to various methods to compute the necessary
        %             calculations               
        %
        %  element_nodes == 2x1 array of node objects that define the nodes
        %                   at the ends of the element

    end
    
    % Public methods go here
    methods (Access = public)
        %% Constructor
        function self = CPBD_Analysis(nnodes,coord,concen,fixity,nele,ends,...
                A,Izz,Iyy,J,Cw,Zzz,Zyy,Ayy,Azz,E,v,Fy,YldSurf,Wt,webdir,...
                beta_ang,w,thermal,truss,anatype)

            self.element_nodes = CreateNodes(self,nnodes,coord);
            self.eleobj = CreateElement (self,A,Ayy,Azz,Iyy,Izz,J,E,v,...
                webdir,ends,nele);

        end

        %% Run the analysis
        function [DEFL, REACT, ELE_FOR, AFLAG] = RunAnalysis(self,nnodes,...
                concen,w,nele,fixity)

            [Pf] = CreateLoadVectors(self,nnodes,concen,w,nele,fixity);

            K = CreateStiffnessMatrix(self,nele,nnodes); 

            [Kff,Kfn,Kfs,Knf,Knn,Kns,Ksf,Ksn,Kss] =...
                ComputeStiffnessSubMatrices(self,K,fixity);

            [self.ELE_FOR] = RecoverElementForces(self,K,nele,fixity,...
                nnodes,concen,w);

            [self.DEFL, self.REACT] =...
                ComputeDisplacementReactions(self,K,fixity,nnodes,concen,...
                w,nele);

            self.AFLAG = CheckKffMatrix(self,Kff);

            [DEFL, REACT, ELE_FOR, AFLAG] = GetMastan2Returns(self,K,...
                nnodes,concen,w,nele,fixity);

        end

        %% Get Mastan2 Returns 
        % Function to return the necessary outputs to Mastan2 to be used by
        % the Mastan2 post processor 
        function [DEFL, REACT, ELE_FOR, AFLAG] = GetMastan2Returns(self,K,...
                nnodes,concen,w,nele,fixity)

            DEFL = self.DEFL
            REACT = self.REACT
            ELE_FOR = self.ELE_FOR
            AFLAG = self.AFLAG;
            ComputeError(self,K,nnodes,concen,w,nele,fixity);

        end
    end
    
    % Private methods go here
    methods (Access = private)

        %% Create Nodes
        % Creates an array of node objects for all of the coordinates given
        % by Mastan in the matrix named "coord"
        function nodes = CreateNodes(self,nnodes,coord)
            nodes(1:nnodes,1) = CPBD_Node([nan,nan,nan]',nan);

            for i = 1:nnodes
                nodes(i) = CPBD_Node(coord(i,:)',i);
            end

        end

        %% Create Element
        % Creates an array of element objects for all the elements given
        % provided by mastan
        function element = CreateElement(self, A, Ayy, Azz, Iyy, Izz, J, E, ...
                v, webdir, ends, nele)
            
            for i = 1:nele
                node1 = self.element_nodes(ends(i,1));
                node2 = self.element_nodes(ends(i,2));
                Nodes = [node1; node2];

                element(i)= CPBD_Element(Nodes, A(i),...
                    Ayy(i), Azz(i), Iyy(i), Izz(i),...
                    J(i), E(i), v(i), webdir(i,:));
            end
        end

        %% Create Stiffness Matrix
        % Creates the entire global stiffness matrix for the structure
        function K = CreateStiffnessMatrix(self,nele,nnodes)

            nDOF = nnodes*6;
            K = zeros(nDOF,nDOF); 

            for i = 1:nele
                K_ele = GetGlobalStiffnessMatrix(self.eleobj(i)); 
                ele_dof = self.eleobj(i).GetElementDOF(); 
                K(ele_dof, ele_dof) = K(ele_dof,ele_dof) + K_ele; 
            end

            K = sparse(K);
        end

        %% Create Load Vectors 
        % Create fixed end force and concentrated force vectors in global
        % coordinates. Also stores Pf which is the concentrated forces
        % minus the fixed end forces acting on the free degrees of freedom
        function [Pf] = CreateLoadVectors(self,nnodes,concen,w,nele,fixity)

            nDOF = nnodes*6;
            P = zeros(nDOF,1);
            concen_trans = (concen)';
        
            for i = 1:nDOF
                P(i,1) = concen_trans(i);
            end

            FeF = zeros(nDOF,1);

            for i=1:nele
                element_dofs = GetElementDOF(self.eleobj(i));
                FeF(element_dofs,1) = FeF(element_dofs,1) +...
                    GetElementFeF(self.eleobj(i),w(i,:).');

            end

            LoadVector = P - FeF;
            [freeDOF, fixedDOF, knownDOF] = ClassifyDOF(self,fixity);
            Pf = LoadVector(freeDOF);

        end

         %% Compute Stiffness SubMatrices
         % Creates the Kff,Kfn, Knf, Knn, Ksf, and Ksn submatricies from ...
         % the global K matrix  
        function [Kff,Kfn,Kfs,Knf,Knn,Kns,Ksf,Ksn,Kss] = ...
                ComputeStiffnessSubMatrices(self,K,fixity)

            [freeDOF, fixedDOF, knownDOF] = ClassifyDOF(self,fixity);

            Kff = K(freeDOF,freeDOF);
            Kfn = K(freeDOF,knownDOF);
            Kfs = K(freeDOF,fixedDOF);
            Knf = Kfn.';
            Knn = K(knownDOF,knownDOF);
            Kns = K(knownDOF,fixedDOF);
            Ksf = Kfs.';
            Ksn = Kns.';
            Kss = K(fixedDOF,fixedDOF);

        end

        %% Classify DOF  
        % Create vectors containing the numbers of the free, displaced, and
        % support degrees of freedom 
        function [freeDOF, fixedDOF, knownDOF] = ClassifyDOF(self,fixity)

            fixity_t = fixity'; 
            freeDOF = find(isnan(fixity_t));
            fixedDOF = find(fixity_t == 0);
            knownDOF = find(fixity_t ~= 0 & isnan(fixity_t)==0);

        end 

        %% Check KFF Matrix 
        % Checks the condition number of the Kff Matrix and reports the
        % number of significant digits lost 
        function flag = CheckKffMatrix(self,Kff)

            k = condest(Kff);
            p = 16;
            s = p - log10(k);

            if s<4
                flag = 0;
            else
                flag = 1;
            end
            disp('The condition number of the Kff matrix is: ')
            disp(k)
            disp('The number of significant digits lost in the analysis is: ')
            disp(log10(k))
        end 

        %% Recover Element Forces 
        % Function that iterates through each element and passes
        % deflections so that element forces are returned in local
        % coordiantes
        function [ELE_FOR] = RecoverElementForces(self,K,nele,fixity,...
                nnodes,concen,w)

            [DEFL, REACT] = ComputeDisplacementReactions(self,K,fixity,...
                nnodes,concen,w,nele);
            element_forces = zeros(nele,12);

            for i = 1:nele 
                element_dof = GetElementDOF(self.eleobj(i));
                DEFL_trans = DEFL.';
                defl_local = DEFL_trans(element_dof);
                element_force = ComputeForces(self.eleobj(i),...
                    defl_local,w(i,:).');
                element_forces(i,1:12) = element_forces(i,1:12) +...
                    element_force.';
            end

            ELE_FOR = element_forces;
        end 

        %% Compute Displacement Reactions
            % Computes the global displacements at free degrees of freedom
            % and formats all displacements for output to MASTAN. Also,
            % computes the global forces at supported and displaced degrees
            % of freedom and formats all the reactions for output to MASTAN

            function [DEFL, REACT] = ComputeDisplacementReactions(self,K,...
                    fixity,nnodes,concen,w,nele)

             [freeDOF, fixedDOF, knownDOF]= ClassifyDOF(self,fixity);
             [Pf] = CreateLoadVectors(self,nnodes,concen,w,nele,fixity);

             [Kff,Kfn,Kfs,Knf,Knn,Kns,Ksf,Ksn,Kss] =...
                 ComputeStiffnessSubMatrices(self,K,fixity);

             fixity_tran = fixity';
             defl_n = fixity_tran(knownDOF);
             defl_f = Kff\(Pf - Kfn*defl_n); 
             defl_s = zeros(length(fixedDOF),1);

             Pn = Knf * defl_f + Knn * defl_n;
             Ps = Ksf * defl_f + Ksn * defl_n;

             defl(freeDOF,1) = defl_f;
             defl(knownDOF,1) = defl_n;
             defl(fixedDOF,1) = defl_s;
             DEFL = zeros(nnodes,6);

             for i = 1:nnodes
                 DEFL(i,1:6) = DEFL(i,1:6) + ...
                     defl(GetNodeDOF(self.element_nodes(i)),1).';
             end

             force = zeros(6*nnodes,1);
             force(freeDOF,1) = Pf;
             force(knownDOF,1) = Pn;
             force(fixedDOF,1) = Ps;

             REACT = zeros(nnodes,6);

             for i = 1:nnodes
                 REACT(i,1:6) = REACT(i,1:6) +...
                     force(GetNodeDOF(self.element_nodes(i)),1).';
             end
        end

        %% Compute Error    
        % Computes the difference between the applied loads Pf and the 
        function ComputeError(self,K,nnodes,concen,w,nele,fixity)

            [Pf] = CreateLoadVectors(self,nnodes,concen,w,nele,fixity);
            [freeDOF, fixedDOF, knownDOF]= ClassifyDOF(self,fixity);
            
            [Kff,Kfn,Kfs,Knf,Knn,Kns,Ksf,Ksn,Kss] = ...
                ComputeStiffnessSubMatrices(self,K,fixity);
            DEFL_transpose = self.DEFL.';
            delF = DEFL_transpose(freeDOF);
            delN = DEFL_transpose(knownDOF);          
            
            Error = Pf - Kff* delF - Kfn* delN;
            disp('Error at free DOFs is:')
            disp(Error)
         end  
    end
end