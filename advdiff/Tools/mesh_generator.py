# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 09:49:18 2022

This module includes all mesh generating functions and mesh quality functions.

@author: Ketil
"""
import logging
import fipy  as fp
import numpy as np


def NSR(vertices):
    """
    Computes the Normalized Shape Ratio for a single cell given coordinates for each vertex.
    This will work for any triangle and rectangular cells, but not general quadrilateral cells.
    The NSR is a float in the range [0,1] where 0 is a very poor value and 1 is an optimal value.

    Parameters
    ----------
    vertices : array ([n_vertices, 2])
               Vertice points (x,y) of a single gridcell. n_vertices = 3 for triangular cells, n_vertices = 4 for rectangular cells.

    Returns
    -------
    rho : float
          The NSR value

    """
    
    if vertices.shape[0] == 3:
        a = np.linalg.norm(vertices[0]-vertices[1])
        b = np.linalg.norm(vertices[1]-vertices[2])
        c = np.linalg.norm(vertices[2]-vertices[0])
        s = (a + b + c) / 2
        num = a * b * c
        den = 2 * np.linalg.norm(np.cross(vertices[0]-vertices[1], vertices[1]-vertices[2]))

        R = num / den                          # Radii of circumcircle of triangle
        r = np.sqrt((s-a) * (s-b) * (s-c) / s) # Radii of incircle of triangle
        
        rho = 2 * r / R  # A perfect equilateral triangle will give rho = 1.0
    
    if vertices.shape[0] == 4:
        c = np.mean(vertices, axis=0)
        dist_x = np.abs(np.min(c[0]-vertices[:,0]))
        dist_y = np.abs(np.min(c[1]-vertices[:,1]))
        
        R = np.max(np.linalg.norm(c-vertices, axis=1)) # Radii of circumcircle of rectangle
        r = np.min([dist_x,dist_y])                    # Radii of incircle of rectangle
        
        rho = np.sqrt(2) * r / R # A perfect square will give rho = 1.0
    
    return rho


def get_NSR(mesh):
    """
    Computes the Normalized Shape Ratio for each cell in a mesh.
    The NSR is a float in the range [0,1] where 0 is a very poor value and 1 is an optimal value.

    Parameters
    ----------
    mesh : fipy.Grid2D 
           A fipy Grid2D object. It can either be a rectilinear or triangular mesh.

    Returns
    -------
    rho : array ([n_cells])
          The NSR value for each cell

    """
    n_cells    = mesh.cellVolumes.shape[0]
    cellIDs    = range(n_cells)
    vertIDs    = mesh._cellVertexIDs.data.T[cellIDs]
    vertCoords = np.array([mesh.vertexCoords.T[vertIDs[ii]] for ii in range(vertIDs.shape[0])])
    rho        = np.array([NSR(vertCoords[ii]) for ii in range(n_cells)])
    
    return rho


def initialize_meshes(grid_type        = 'equi',
                      grid_type_out    = 'equi',
                      Lx               = 10000.,
                      Ly               = 10000.,
                      nx               = 100,
                      ny               = 100,
                      minvol           = 10000,
                      maxvol           = 100,
                      reduction_factor = 0.96,
                      power            = 1.0,
                      origin           = np.array([0.,0.]), 
                      verbose          = True):
    """
    Initializes the computational and output meshes for use in FiPy. 

    Parameters
    ----------
    grid_type         : str 
                        Type of computational grid. 
    grid_type_out     : str
                        Type of output grid.
    Lx/Ly             : float
                        Length of the domain in the x/y direction. 
    nx/ny             : int
                        Number of discretizations in x/y direction in a rectilinear mesh. 
    maxvol            : float
                        The maximum volume for the largest gridcell in a triangular mesh.
    minvol            : float
                        The minimum volume for the smallest gridcell in a triangular mesh.
    reduction_factor  : float
                        How much strong the exponential reduction of gridcell size is when using expnential grid.
    power             : float
                        A float greater or equal to zero. Only relevant for triangular meshes. Higher power means more gridcells. By default 1.0.
                        Power = 0.0 --> cell volumes are constant = maxvol.
                        Power = 0.5 --> cell volumes increase proportionally to the square root of the distance from the center from minvol to maxvol.
                        Power = 1.0 --> cell volumes increase linearily by the distance from the center from minvol to maxvol.
                        Power = 2.0 --> cell volumes increase proportionally to the square of the distance from the center from minvol to maxvol.
    origin            : array ([x,y])
                        Where to set the origin (the middle of the mesh) in the global domain.
    verbose           : bool
                        If True, print mesh information.

    Returns
    -------
    meshdict        : dict
                      Dictionary containing two fipy Grid2D objects and two boolean variables
    mesh            : fipy.Grid2D
                      computational mesh 
    mesh_out        : fipy.Grid2D 
                      output mesh        
    is_onemesh      : bool
                      True if mesh == mesh_out else False
    is_unstructured : bool
                      True if mesh is unstructured

    """

    # If the comp mesh and the out mesh are the same, then we only generate one mesh.
    # Therefore there is no need for interpolation. The solver will be aware of this if is_onemesh == True.
    if grid_type_out == grid_type:
        is_onemesh  = True
        mesh     = create_mesh(Lx=Lx, Ly=Ly, nx=nx, ny=ny, maxvol=maxvol, minvol=minvol, 
                               reduction_factor=reduction_factor, power=power, grid_type=grid_type, origin=origin)
        mesh_out = mesh
    else: # Else the solver needs to interpolate between the two different meshes...
        is_onemesh  = False
        mesh     = create_mesh(Lx=Lx, Ly=Ly, nx=nx, ny=ny, maxvol=maxvol, minvol=minvol, 
                               reduction_factor=reduction_factor, power=power, grid_type=grid_type, origin=origin)
        mesh_out = create_mesh(Lx=Lx, Ly=Ly, nx=nx, ny=ny, maxvol=maxvol, minvol=minvol, 
                               reduction_factor=reduction_factor, power=power, grid_type=grid_type_out, origin=origin)
        
    is_unstructured = True if grid_type == 'tri' else False 
        
    if verbose:
        vols     = mesh.cellVolumes
        vols_out = mesh_out.cellVolumes
        NSR_score     = get_NSR(mesh=mesh)
        NSR_out_score = get_NSR(mesh=mesh_out)
        logging.info('--- Computational mesh info ---')
        logging.info(f'Number of Cells:       {vols.shape[0]}')
        logging.info(f'Volume  0.1% Q (min):  {np.quantile(vols, 0.001):.2f} ({np.min(vols):.2f})')
        logging.info(f'Volume 50.0% Q (mean): {np.quantile(vols, 0.500):.2f} ({np.mean(vols):.2f})')
        logging.info(f'Volume 99.9% Q (max):  {np.quantile(vols, 0.999):.2f} ({np.max(vols):.2f})')
        logging.info(f'NSR min, mean, max:    {np.min(NSR_score):.2f}, {np.mean(NSR_score):.2f}, {np.max(NSR_score):.2f}')
        logging.info('--- Output mesh info ---')
        logging.info(f'Number of Cells:       {vols_out.shape[0]}')
        logging.info(f'Volume  0.1% Q (min):  {np.quantile(vols_out, 0.001):.2f} ({np.min(vols_out):.2f})')
        logging.info(f'Volume 50.0% Q (mean): {np.quantile(vols_out, 0.500):.2f} ({np.mean(vols_out):.2f})')
        logging.info(f'Volume 99.9% Q (max):  {np.quantile(vols_out, 0.999):.2f} ({np.max(vols_out):.2f})')
        logging.info(f'NSR min, mean, max:    {np.min(NSR_out_score):.2f}, {np.mean(NSR_out_score):.2f}, {np.max(NSR_out_score):.2f}')
        logging.info('Meshes generated...\n')
    
    meshdict = {
        'mesh'            : mesh,
        'mesh_out'        : mesh_out,
        'mesh_type'       : grid_type,
        'mesh_out_type'   : grid_type_out,
        'is_onemesh'      : is_onemesh,
        'is_unstructured' : is_unstructured,
        'Lx'              : Lx,
        'Ly'              : Ly,
        'nx'              : nx,
        'ny'              : ny,
        'maxvol'          : maxvol,
        'minvol'          : minvol,
        'reduction_factor': reduction_factor,
        'power'           : power,
        'origin'          : origin,
        }
    
    return meshdict


def tri_mesh(Lx, Ly, maxvol=25000, minvol=250, origin=(0,0), iterations=3, power=1.0):
    """
    Custom triangle mesh using Gmsh.
    Generates a fipy.Grid2D object containing a triangle mesh constructed using GMSH. 
    The mesh has high density in the middle and low density along the boundaries.

    See https://www.ctcms.nist.gov/fipy/download/fipy-3.1.3.pdf on page 204 and below for more info.

    Parameters
    ----------
    Lx : float
         Length of domain along x-axis.
    Ly : float
         Length of domain along y-axis.
    maxvol : float
             Maximum volume of largest  gridcell in triangle mesh (Approximate)
    minvol : float
             Minimum volume of smallest gridcell in triangle mesh (Approximate)
    origin : array ([x,y])
             Where to set the origin (the middle of the mesh) of domain.
    iterations : int
                 Number of refinement iterations to spend refining the mesh. (3-4 is enough.)
    power : float
            A float greater or equal to zero. Higher power means more gridcells. By default 1.0.
            Power = 0.0 --> cell volumes are constant = maxvol.
            Power = 0.5 --> cell volumes increase proportionally to the square root of the distance from the center from minvol to maxvol.
            Power = 1.0 --> cell volumes increase linearily by the distance from the center from minvol to maxvol.
            Power = 2.0 --> cell volumes increase proportionally to the square of the distance from the center from minvol to maxvol.

    Returns
    -------
    mesh : fipy.Grid2D
           A fipy.Grid2D object containing a refined triangle mesh.

    """

    # Define GMSH geometry script
    geo = '''

    // A mesh consisting of a square


    // define the corners of the square
    
    
    Point(1) = {%f, %f, 0, 1};
    
    Point(2) = {%f, %f, 0, 1};
    
    Point(3) = {%f, %f, 0, 1};
    
    Point(4) = {%f, %f, 0, 1};
    
    
    // define the square
    
    
    Line(1) = {1, 2};
    
    Line(2) = {2, 3};
    
    Line(3) = {3, 4};
    
    Line(4) = {4, 1};
    
    
    // define the boundary
    
    
    Line Loop(1) = {1, 2, 3, 4};
    
    
    // define the domain
    
    
    Plane Surface(1) = {1};
    
    Mesh.MeshSizeExtendFromBoundary = 0;
    Mesh.MeshSizeFromPoints = 0;
    Mesh.MeshSizeFromCurvature = 0;
    
    '''
    
    geo = geo % (Lx/2, Ly/2, -Lx/2, Ly/2, -Lx/2, -Ly/2, Lx/2, -Ly/2) # Scale unit rectangle
    

    def sidelength_func(x,y):              # Sidelength function. Minima at origin.
        a = np.sqrt(4*maxvol/np.sqrt(3))   # Sidelenghts of an equilateral triangle given some approximate maximum area
        b = np.sqrt(4*minvol/np.sqrt(3))   # Sidelengths of an equilateral triangle given some approximate minimum area
        return (a-b)*(np.sqrt(((x)**2+(y)**2)/((Lx/2)**2+(Ly/2)**2)))**power+b

    def sidelength_func_OLD(x,y):          # Sidelength function. Minima at origin.
        a = 2*((2+np.sqrt(2))/3)*maxvol    # Sidelenghts of an equilateral triangle given some approximate maximum area (Incorrect)
        b = 2*((2+np.sqrt(2))/3)*minvol    # Sidelengths of an equilateral triangle given some approximate minimum area (Incorrect)
        d = (np.sqrt(a)-np.sqrt(b))**2 - a # Correction term
        return np.sqrt(a+d)*np.sqrt(((x)**2+(y)**2)/((Lx/2)**2+(Ly/2)**2))+np.sqrt(b)

    bkg = None

    for refine in range(iterations): # Converge onto a refined mesh over some iterations

        mesh = fp.Gmsh2D(geo, background=bkg) 
        x, y = mesh.cellCenters 

        # Average triangle sidelength given by the value of the cone surface
        bkg = fp.CellVariable(mesh=mesh, value=sidelength_func(x,y))
        #bkg = fp.CellVariable(mesh=mesh, value=sidelength_func_OLD(x,y))

    mesh = fp.Gmsh2D(geo, background=bkg)     # Generate mesh using final refinement bkg
    mesh = mesh + ((origin[0],),(origin[1],)) # Translate mesh to relevant origin
    return mesh


def create_mesh(Lx=10000, Ly=10000, nx=100, ny=100, maxvol=25000, minvol=250, reduction_factor=0.96, power=1.0,
                grid_type='equi', origin=(0,0)):
    '''
    Constructs a rectangular mesh with the origin at the center at location given by origin=[x,y].
    The grid can be constructed in one of four ways:

    *   grid_type='equi':
        The distance between gridpoints is equidistant.
    *   grid_type='segmented':
        The distance between gridpoints is hard coded, for which the distance
        between gridpoints are larger near the edges and small around the origin.
    *   grid_type='exp':
        The distance between gridpoints reduces exponentially from the edges to the origin.
        The coefficient of change is determined by the reduce factor.
    *   grid_type='tri':
        The domain is divided into a mesh of triangles for which the density increases towards the origin.
        The minimum volume of the smallest cells is set by minvol, while the maximum of the largest cells volume is set by maxvol.
    
    Parameters
    ----------
    grid_type : str
                Type of mesh to generate (Allowable: 'equi', 'exp', 'segmented', 'tri')
    Lx        : float
                Length of domain along x-axis.
    Ly        : float
                Length of domain along y-axis.
    nx        : int
                Number of meshpoints along x-axis (Only relevant for 'equi', 'segmented' and 'exp').
    ny        : int
                Number of meshpoints along y-axis (Only relevant for 'equi', 'segmented' and 'exp').
    origin    : array ([x,y])
                Where to set the origin (the middle of the mesh) of domain.
    maxvol    : float
                Maximum volume of largest  gridcell in triangle mesh (Approximate) (Only relevant for 'tri')
    minvol    : float
                Minimum volume of smallest gridcell in triangle mesh (Approximate) (Only relevant for 'tri')
    reduction_factor : float
                       How much strong the exponential reduction of gridcell size is when using expnential grid.
    power     : float
                A float greater or equal to zero. Higher power means more gridcells. By default 1.0.
                Power = 0.0 --> cell volumes are constant = maxvol.
                Power = 0.5 --> cell volumes increase proportionally to the square root of the distance from the center from minvol to maxvol.
                Power = 1.0 --> cell volumes increase linearily by the distance from the center from minvol to maxvol.
                Power = 2.0 --> cell volumes increase proportionally to the square of the distance from the center from minvol to maxvol.

    Returns
    -------
    fipy.Grid2D object.

    '''
    
    if grid_type == 'tri':
        return tri_mesh(Lx=Lx, Ly=Ly, maxvol=maxvol, minvol=minvol, origin=origin, power=power)
        
    else:
        def dxvec_equidistant(nx):
            dxvec = np.array([1.0/nx]*nx)
            return dxvec
    
        def dxvec_nonuniform_segments(nx):
            L1, L2, L3 = 100, 10, 100
            n1, n3 = int(0.1*nx), int(0.1*nx)
            n2 = nx - n1 - n3
            dx1, dx2, dx3 = L1/n1, L2/n2, L3/n3
            dxvec = np.array([dx1]*n1 + [dx2]*n2 + [dx3]*n3)
            return dxvec
    
        def dxvec_nonuniform_exp(nx, reduction_factor):
            reduce_factor = reduction_factor # Should be 0.96 
            ndwn=int(np.floor(nx/2))
            nup=int(np.ceil(nx/2))            
            d1 = np.array([reduce_factor**i for i in range(ndwn)])
            d2=np.array([reduce_factor**i for i in range(nup)])
            dxvec=np.hstack((d1,np.flipud(d2)))
            dxvec=dxvec/np.sum(dxvec)
            return dxvec
    
        if grid_type == 'equi':
            dx = dxvec_equidistant(nx)*Lx
            dy = dxvec_equidistant(ny)*Ly
        
        elif grid_type == 'exp':
            dx = dxvec_nonuniform_exp(nx, reduction_factor)*Lx 
            dy = dxvec_nonuniform_exp(ny, reduction_factor)*Ly 
    
        elif grid_type == 'segmented':
            dx = dxvec_nonuniform_segments(nx)*Lx 
            dy = dxvec_nonuniform_segments(ny)*Ly
        
        
        return fp.Grid2D(dx=dx, dy=dy) - ((Lx/2,),(Ly/2,)) + ((origin[0],),(origin[1],))

  