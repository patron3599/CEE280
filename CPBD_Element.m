classdef CPBD_Element < handle

% This class contains the necesary properties and methods to define the
% framing (i.e beams and columns) parameters of the structural model. When
% this class is called, an instance of this class will be created for each
% element in the structure. This instance will contain all relevant
% properties such as its length, nodes, stiffness matrix, transformation
% matrix, and its degrees of freedom. This instance can then be called from
% other classes to run the analysis.
    
    % Private for the element class 
    properties (Access = private) 
        
        length; 
        Elastic_Stiffness_Matrix; 
        Transformation_Matrix; 
        element_dof;
        element_nodes;

        %% Property Definitions 
        %  Note: Variables defined in ud_3d1el.m are not repeated here
        %
        %  length  == Element's length
        %
        %  Elastic_Stiffness_Matrix == Element's global stiffness matrix
        %
        %  Transformation_Matrix == Element's global transformation matrix 
        %
        %  element_dof == 12x1 vector containing the numbers of degree of 
        %                 freedom for the element in the order of: delta_x,
        %                 delta_y, delta_z, theta_x, theta_y, theta_z. This
        %                 also applies to the other end of the element
        %
        %  element_nodes == 2x1 array of node objects that define the nodes
        %                   at the ends of the element

    end
    
    % Public methods go here
    methods (Access = public)
        %% Constructor
        function self = CPBD_Element(element_nodes, A, Ayy, Azz, Iyy, Izz, ...
                J, E, v, webdir)

                self.element_nodes = element_nodes;
                self.length = ComputeLength(self, element_nodes);

                self.Transformation_Matrix = ...
                    ComputeTransformationMatrix(self, element_nodes,...
                    webdir);

                self.Elastic_Stiffness_Matrix = ...
                    ComputeElasticStiffnessMatrix(self, A, Ayy, Azz,...
                    Iyy, Izz, J, E, v);   

                self.element_dof=RetrieveDOF(self);
        end

    %% Get Global Stiffness Matrix
    % Function to query global stiffness matrix from outside the element
    % class
    function Elastic_Stiffness_Matrix = GetGlobalStiffnessMatrix(self)
        Elastic_Stiffness_Matrix = self.Elastic_Stiffness_Matrix;
    end

    %% Get Element DOF's
    % Funtion that assins 12x1 vector property of element class to output 
    function element_dof = GetElementDOF(self)
        element_dof = self.element_dof;
    end

    %% Compute Forces 
    % Computes the local element forces based on the element's nodal
    % deflections 
    function element_forces = ComputeForces(self, defl_local, w) 

        K_Local = self.Transformation_Matrix * ...
            self.Elastic_Stiffness_Matrix * self.Transformation_Matrix';

        element_forces = K_Local* ...
            self.Transformation_Matrix * defl_local + ...
            ComputeFixedEndForces(self,w);

    end

    %% Get the global FEF's
    % Function to query global stiffness matrix from outside the element
    % class
        function FeF_Global = GetElementFeF(self,w)

            FeF_local = ComputeFixedEndForces(self,w);
            FeF_Global = self.Transformation_Matrix' * FeF_local;
           
        end
        end
  
    
    % Private methods go here
    methods (Access = private)
        %% Compute the element's length
        % Computes the length of a single element 
        function length = ComputeLength(self, element_nodes)

            node1 = element_nodes(1,:);
            node2 = element_nodes(2,:);
            node1_coord = GetNodeCoord(node1);
            node2_coord = GetNodeCoord(node2);

            length = norm(node2_coord-node1_coord);
 
        end
        
        %% Compute the element's geometric transformation matrix
        % This function computes the geometric transformation matrix which
        % allows for the properties of the element to be converted from
        % local to global coordinates and vice versa.
        function Transformation_Matrix = ComputeTransformationMatrix(self,...
                element_nodes, webdir)

            node1 = element_nodes(1,:);
            node2 = element_nodes(2,:);

            node1_coord = GetNodeCoord(node1);
            node2_coord = GetNodeCoord(node2);

            x_prime = (node2_coord - node1_coord);
            x_prime = x_prime/norm(x_prime);
            y_prime = webdir';
            y_prime = y_prime/norm(y_prime);
            z_prime = cross(x_prime,y_prime); 
            z_prime = z_prime/norm(z_prime);
            gamma_prime = [x_prime, y_prime, z_prime]';

            Transformation_Matrix = blkdiag(gamma_prime,gamma_prime,...
                gamma_prime,gamma_prime);

        end
        
        %% Compute the element's Elastic Stiffness Matrix 
        % This function computes the element's stiffness matrix in local
        % and global coordinates
        function Elastic_Stiffness_Matrix = ...
                ComputeElasticStiffnessMatrix(self, A, Ayy, Azz, Iyy, Izz, ...
                J, E, v)

            % Compute missing parameters
                G = E/(2*(1+v)); % Modulus of Rigidity 

                phi_z = (12*E*Iyy)/(G*Azz*self.length...
                    ^(2)); % Z-Direction Shear deformation parameter

                phi_y = (12*E*Izz)/(G*Ayy*self.length...
                    ^(2));% Y-Direction Shear deformation parameter
                
            % Construct Local Element Stiffness Matrix 
            k_11 = (E*A)/self.length;
            k_17 = -k_11;
            k_22 = (12*E*Izz)/((self.length^(3))*(1+phi_y));
            k_26 = (6*E*Izz)/(self.length^(2)*(1+phi_y));
            k_28 = -k_22;
            k_212 = k_26;
            k_33 = (12*E*Iyy)/((self.length^(3))*(1+phi_z));
            k_35 = (-6*E*Iyy)/(self.length^(2)*(1+phi_z));
            k_39 = (-12*E*Iyy)/((self.length^(3))*(1+phi_z));
            k_311 = (-6*E*Iyy)/(self.length^(2)*(1+phi_z));
            k_44 = (G*J)/self.length;
            k_410 = -k_44;
            k_53 = k_35;
            k_55 = ((4+phi_z)*E*Iyy)/(self.length*(1+phi_z));
            k_59 = -k_53;
            k_511 = ((2-phi_z)*E*Iyy)/(self.length*(1+phi_z));
            k_62 = k_26;
            k_66 = ((4+phi_y)*E*Izz)/(self.length*(1+phi_y));
            k_68 = -k_62;
            k_612 = ((2-phi_y)*E*Izz)/(self.length*(1+phi_y));
            k_71 = -k_11;
            k_77 = -k_17;
            k_82 = -k_22;
            k_86 = -k_26;
            k_88 = -k_28;
            k_812 = -k_212;
            k_93 = -k_33;
            k_95 = -k_35;
            k_99 = -k_39;
            k_911 = -k_311;
            k_104 = -k_44;
            k_1010 = -k_410;
            k_113 = -k_53;
            k_115 = k_511;
            k_119 = k_59;
            k_1111 = k_55;
            k_122 = k_62;
            k_126 = k_612;
            k_128 = k_68;
            k_1212 = k_66;
            
            k_ele_3d = [  k_11    0     0     0     0     0    k_17   0     0     0      0       0;
                           0     k_22   0     0     0    k_26   0    k_28   0     0      0      k_212;
                           0      0    k_33   0    k_35   0     0     0    k_39   0     k_311    0;
                           0      0     0    k_44   0     0     0     0     0    k_410   0       0;
                           0      0    k_53   0    k_55   0     0     0    k_59   0     k_511    0;
                           0     k_62   0     0     0    k_66   0    k_68   0     0      0      k_612;
                          k_71    0     0     0     0     0    k_77   0     0     0      0       0;
                           0     k_82   0     0     0    k_86   0    k_88   0     0      0      k_812;
                           0      0    k_93   0    k_95   0     0     0    k_99   0     k_911    0;
                           0      0     0    k_104  0     0     0     0     0    k_1010  0       0;
                           0      0    k_113  0    k_115  0     0     0    k_119   0    k_1111   0;
                           0     k_122  0     0     0    k_126  0    k_128  0      0     0      k_1212  ];

            % Compute Global Elastic Stiffness Matrix
            Elastic_Stiffness_Matrix = self.Transformation_Matrix'*...
                k_ele_3d*self.Transformation_Matrix;
            
        end
        
         %% Compute Fixed End Forces  
            % This function is used to calculate the Fixed End Forces of an
            % element
            function FeF_local = ComputeFixedEndForces(self,w)

                w_x = w(1,:);
                w_y = w(2,:);
                w_z = w(3,:);
         
                FeF_1 = (-w_x*self.length)/2;
                FeF_7 = FeF_1; 
                FeF_2 = (-w_y*self.length)/2;
                FeF_8 = FeF_2;
                FeF_3 = (-w_z*self.length)/2;
                FeF_9 = FeF_3;
                FeF_4 = [0];
                FeF_10 = FeF_4;
                FeF_5 = (w_z*self.length^(2))/12;
                FeF_11 = -FeF_5;
                FeF_6 = -(w_y*self.length^(2))/12;
                FeF_12 = -FeF_6;

        FeF_local = [FeF_1;FeF_2;FeF_3;FeF_4;FeF_5;FeF_6;FeF_7;...
            FeF_8;FeF_9;FeF_10;FeF_11;FeF_12];

            end

    %% Retrieve DOF's
    % This function is used to assign degress of freedom to an element....
    % It is then placed in our constructor of this class to then be ...
    % called for other calculations.   
    function element_dof = RetrieveDOF(self)

        node1 = self.element_nodes(1,:);
        node2 = self.element_nodes(2,:);

        element_node1_dof = GetNodeDOF(node1);
        element_node2_dof = GetNodeDOF(node2);

        element_dof = [element_node1_dof; element_node2_dof]; 

    end

    end
    end
