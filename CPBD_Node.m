classdef CPBD_Node < handle

% This class contains the necesary properties and methods to define the
% nodal parameters of the structural model. Each instance of the Node Class 
% will contain all relevant nodal properties (i.e coordinates and degrees 
% of freedom) and methods that can be called from other classes that are 
% required to run the analysis.

    
    % Private properties go here
    properties (Access = private)

        node_coord; 
        node_number;

        %% Property Definitions 
        %  Note: Variables defined in ud_3d1el.m are not repeated here
        %
        %  node_coord  == 3x1 vector containing the x, y, and z coordinates
        %                of the node 
        %
        %  node_number == User defined value corresponding to its location
        %                on the structure 
   
    end
    
    % Public methods go here
    methods (Access = public)
        %% Constructor
        function self = CPBD_Node(node_coord, node_number)
            self.node_coord = node_coord;
            self.node_number = node_number;

        end
        
        %% Get Node Coordinates
        %  Returns node's coordiantes for access outside the node class
        function node_coord = GetNodeCoord(self)
            node_coord = self.node_coord;
        end

        %% Get Node DOFs
        %  Returns the node's 6 DOF's for access outside the node class 
        function node_dof = GetNodeDOF(self)
            node_dof = AssignDOF(self);
        end
    end
    
    % Private methods go here
    methods (Access = private)
        %% Assign DOF Numbers 
        %   Assigns the DOF numbers corresponding to each node number 
        function node_dof = AssignDOF(self)
            node_dof = 6*(self.node_number-1) + (1:6)'; 
        end
    end
end
