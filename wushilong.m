function exam8_1
% This program is a calculation program in Chapter 8 of "Structural Analysis Finite Element Method and MATLAB Design". 
%It uses plane beam elements to calculate the free vibration characteristics of two-hinge parabolic arches.
%      Input parameters: none
%      Output: The first 3 vibration frequencies and their corresponding mode shapes

% define global variables
%      gNode ------ Node coordinates
%      gElement --- unit definition
%      gMaterial -- Material properties
%      gBC1 ------- The first kind of constraints
%      gK --------- Global stiffness matrix
%      gDelta ----- Overall Node Coordinates

    PlaneFrameModel ;      % Defining the Finite Element Model
    SolveModel ;           % Solve the Finite Element Model
    DisplayResults ;       % Display calculation results
return ;

function PlaneFrameModel
%  Defining a finite element model of a planar rod system
%  Input parameters:
%      none
%  Return value:
%      none
%  instruction:
%      This function defines the finite element model data of the planar rod system:
%        gNode ------- Node definition
%        gElement ---- unit definition
%        gMaterial --- Material definitions, including elastic modulus, beam cross-sectional area, and beam bending moment of inertia
%        gBC --------- Restrictions

    global gNode gElement gMaterial gBC1

    % geometric features of a given parabolic arch
    L = 60 ;               %  Calculate span  
    f = 7.5 ;              %  Calculate the sag
    
    n = 100 ;              %  number of units
    x = -L/2:L/n:L/2 ;     %  the x coordinate of the node
    a = f/L^2*4 ;
    y = - a * x.^2 ;       %  the y coordinate of the node

    % Node coordinates
    gNode = [x'  y'] ;
    
    % unit definition
    gElement = zeros( n, 3 ) ;
    for i=1:n
        gElement( i, : ) = [ i, i+1, 1 ] ;
    end
    
    % Material properties
    %Elastic Modulus  bending moment of inertia  Cross-sectional area  density
    gMaterial = [2.06e11,  0.03622,   0.0815,  1435.2/0.0815];   

    % The first kind of constraints
    % node number /degrees of freedom  /Constraint value
    gBC1 = [ 1,        1,        0.0
             1,        2,        0.0
             n+1,      1,        0.0
             n+1,      2,        0.0] ;
return

function SolveModel
%  Solve the Finite Element Model
%  Input parameters:
%     none
%  Return value:
%     none


    global gNode gElement gMaterial gBC1 gK gM gEigValue gEigVector

    % step1. Define global stiffness matrix and nodal force vectors
    [node_number,dummy] = size( gNode ) ;
    gK = sparse( node_number * 3, node_number * 3 ) ;
    gM = sparse( node_number * 3, node_number * 3 ) ;

    % step2. Calculate element stiffness and mass matrices and integrate into global stiffness and mass matrices
    [element_number,dummy] = size( gElement ) ;
    for ie=1:1:element_number
        k = StiffnessMatrix( ie ) ;
        m = MassMatrix( ie ) ; 
        AssembleGlobalMatrix( ie, k, m ) ;
    end

    % step3. To deal with the first type of constraints, modify the stiffness and mass matrices (using line-by-line method)
    [bc1_number,dummy] = size( gBC1 ) ;
    w2max = max( diag(gK)./diag(gM) ) ;
    for ibc=1:1:bc1_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*3 + d ;
        gK(:,m) = zeros( node_number*3, 1 ) ;
        gK(m,:) = zeros( 1, node_number*3 ) ;
        gK(m,m) = 1;
        gM(:,m) = zeros( node_number*3, 1 ) ;
        gM(m,:) = zeros( 1, node_number*3 ) ;
        gM(m,m) = gK(m,m)/w2max/1e10 ;
    end
    
    % step4. Solve the eigenvalue problem
    % step4.1In order to make the stiffness and mass matrices symmetric (round-off errors may be introduced in the calculation)
    for i=1:node_number*3
        for j=i:node_number*3
            gK(j,i) = gK(i,j) ;
            gM(j,i) = gM(i,j) ;
        end
    end
    
    % step4.2 Calculate the first 6 order eigenvalues and eigenvectors
    [gEigVector, gEigValue] = eigs(gK, gM, 3, 'SM' ) ; 
    
    % step4.3 Modify constrained degrees of freedom in eigenvectors
    for ibc=1:1:bc1_number
        n = gBC1(ibc, 1 ) ;
        d = gBC1(ibc, 2 ) ;
        m = (n-1)*3 + d ;
        gEigVector(m,:) = gBC1(ibc,3) ;
    end
return

function k = StiffnessMatrix( ie )
%  Calculate the element stiffness matrix
%  Input parameters:
%     ie -------  unit number
%  Return value:
%     k  ----  Stiffness Matrix in Global Coordinate System
    global gNode gElement gMaterial
    k = zeros( 6, 6 ) ;
    E = gMaterial( gElement(ie, 3), 1 ) ;
    I = gMaterial( gElement(ie, 3), 2 ) ;
    A = gMaterial( gElement(ie, 3), 3 ) ;
    xi = gNode( gElement( ie, 1 ), 1 ) ;
    yi = gNode( gElement( ie, 1 ), 2 ) ;
    xj = gNode( gElement( ie, 2 ), 1 ) ;
    yj = gNode( gElement( ie, 2 ), 2 ) ;
    L = ( (xj-xi)^2 + (yj-yi)^2 )^(1/2) ;
    k = [  E*A/L           0          0 -E*A/L           0          0
               0  12*E*I/L^3  6*E*I/L^2      0 -12*E*I/L^3  6*E*I/L^2
               0   6*E*I/L^2    4*E*I/L      0  -6*E*I/L^2    2*E*I/L
          -E*A/L           0          0  E*A/L           0          0
               0 -12*E*I/L^3 -6*E*I/L^2      0  12*E*I/L^3 -6*E*I/L^2
               0   6*E*I/L^2    2*E*I/L      0  -6*E*I/L^2    4*E*I/L] ;
    T = TransformMatrix( ie ) ;
    k = T*k*transpose(T) ;
return

function m = MassMatrix( ie )
%  Calculate the element mass matrix
%  Input parameters:
%     ie -------  unit number
%  Return value:
%     m  ----  Mass matrix in global coordinate system
    global gNode gElement gMaterial
    m = zeros( 6, 6 ) ;
    E = gMaterial( gElement(ie, 3), 1 ) ;
    A = gMaterial( gElement(ie, 3), 3 ) ;
    ro = gMaterial( gElement(ie, 3 ), 4 ) ;
    xi = gNode( gElement( ie, 1 ), 1 ) ;
    yi = gNode( gElement( ie, 1 ), 2 ) ;
    xj = gNode( gElement( ie, 2 ), 1 ) ;
    yj = gNode( gElement( ie, 2 ), 2 ) ;
    L = ( (xj-xi)^2 + (yj-yi)^2 )^(1/2) ;
    m = ro*A*L/420*[140      0      0   70      0      0
                      0    156   22*L    0     54   -13*L
                      0   22*L  4*L^2    0   13*L  -3*L^2
                     70      0      0  140      0       0 
                      0     54   13*L    0    156   -22*L
                      0  -13*L -3*L^2    0  -22*L  4*L^2 ] ;
    T = TransformMatrix( ie ) ;
    m = T*m*transpose(T) ;
return

function AssembleGlobalMatrix( ie, ke, me )
%  Integrate element stiffness and mass matrices into the global stiffness matrix
%  Input parameters:
%      ie  --- unit number
%      ke  --- Element Stiffness Matrix
%      me  --- Element Mass Matrix
%  return value:
%     none
    global gElement gK gM
    for i=1:1:2
        for j=1:1:2
            for p=1:1:3
                for q =1:1:3
                    m = (i-1)*3+p ;
                    n = (j-1)*3+q ;
                    M = (gElement(ie,i)-1)*3+p ;
                    N = (gElement(ie,j)-1)*3+q ;
                    gK(M,N) = gK(M,N) + ke(m,n) ;
                    gM(M,N) = gM(M,N) + me(m,n) ;
                end
            end
        end
    end
return

function T = TransformMatrix( ie )
%  Calculate the coordinate transformation matrix of the unit (local coordinates -> global coordinates)
%  Input parameters
%      ie  ----- node number
%  return value
%      T ------- Coordinate transformation matrix from local coordinates to global coordinates
    global gElement gNode
    xi = gNode( gElement( ie, 1 ), 1 ) ;
    yi = gNode( gElement( ie, 1 ), 2 ) ;
    xj = gNode( gElement( ie, 2 ), 1 ) ;
    yj = gNode( gElement( ie, 2 ), 2 ) ;
    L = sqrt( (xj-xi)^2 + (yj-yi)^2 ) ;
    c = (xj-xi)/L ;
    s = (yj-yi)/L ;
    T=[ c  -s   0   0   0   0
        s   c   0   0   0   0
        0   0   1   0   0   0
        0   0   0   c  -s   0
        0   0   0   s   c   0
        0   0   0   0   0   1] ;
return

function DisplayResults
%  Display calculation results
%  Input parameters:
%      none
%  Return value:
%     none

    global gNode gElement gMaterial gBC1 gEigValue gEigVector

    fre_number = length(diag(gEigValue)) ;
    
    % print eigenvectors (mode shapes)
    fprintf( '\n\n  Table 1 Eigenvectors (Mode Shapes)  \n' ) ;
    for i=1:fre_number
        fprintf( '----------------') ;
    end
    fprintf( '\n' ) ;
    for i=1:fre_number
        fprintf( '  %6d        ', i ) ;
    end
    fprintf( '\n' ) ;
    for i=1:fre_number
        fprintf( '----------------') ;
    end
    fprintf( '\n' ) ;
    [dof,dummy]=size(gEigVector) ;
    for i=1:dof
        for j=fre_number:-1:1
            fprintf( '%15.7e ', gEigVector(i,j) ) ;
        end
        fprintf( '\n' ) ;
    end
    for i=1:fre_number
        fprintf( '----------------') ;
    end
    fprintf( '\n' ) ;

    % print eigenvalues
    fprintf( '\n\n\n\n Table 2 List of eigenvalues (frequency)  \n' ) ;
    fprintf( '-----------------------------------------------------------------------\n') ;
    fprintf( '   Order    Eigenvalues      Frequency (Hz)   Circular frequency (Hz) \n' ) ;
    fprintf( '-----------------------------------------------------------------------\n') ;
    for i=fre_number:-1:1
        fprintf( '%6d   %15.7e   %15.7e   %15.7e\n', fre_number-i+1, ...
            gEigValue(i,i), sqrt(gEigValue(i,i))/2/pi, sqrt(gEigValue(i,i)) ) ;
    end
    fprintf( '----------------------------------------------------------------------\n') ;
    
    % Plot Mode Shapes
    for j=fre_number:-1:1
        figure ;
        x = gNode(:,1) ;
        y = gNode(:,2) ;
        dx = gEigVector(1:3:length(x)*3, j ) ;
        dy = gEigVector(2:3:length(x)*3, j ) ;
        factor = max( [max(abs(x))/max(abs(dx)), max(abs(y))/max(abs(dy))] )* 0.05; 
        plot(x,y,'-', x+factor*dx, y+factor*dy, ':') ;
        title( sprintf( '%d order frequency: %.3f Hz', fre_number-j+1, sqrt(gEigValue(j,j))/2/pi ) ) ;
        axis equal;
        axis off ;
    end
return
