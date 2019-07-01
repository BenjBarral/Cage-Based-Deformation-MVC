# An interactive GUI to animate character meshes using cage based deformation.

![](media/giphy.gif)

An application that allows users to deform meshes in a smooth and realistic way in real-time. Based on cage-based deformations : the user moves a coarse version of the mesh (a "cage") and the mesh deforms by using vertex-to-mesh interpolation (mean value coordinates).
You can either use existing (hand made) cage meshes, or the program can generate a coarse cage automatically.

Implementetation of two papers:
- Automatic coarse cage generation : 
Xian et al. Automatic Generation of Coarse Bounding Cages from Dense Meshes ([paper](
http://www.cad.zju.edu.cn/home/hwlin/pdf_files/Automatic-generation-of-coarse-bounding-cages-from-dense-meshes.pdf))
- Mean value coordinates (for vertex-cage interpolation) : 
Ju et al. Mean Value Coordinates for Closed Triangular Meshes ([paper](https://www.cse.wustl.edu/~taoju/research/meanvalue.pdf))


## Usage
### Dependencies
You need the [LibIGL](https://github.com/libigl/libigl) library.
Put the libigl source code folder two levels above this repository.

### Building and running the code
Build using CMake and the CMakeLists.txt provided.
Run the ```3DA_project_CageBasedDef``` target.
Either use the provided mesh and cage files (or your own), or generate the cage with code using the cage generation algorithm: set the boolean ```compute_automatic_cage``` accordingly.
Change the path the mesh you want to use in the ```mesh_file_name``` variable.

### Using the interface
Deform the cage by clicking and dragging the handle points (see the video in intro for an example of user action) 
- Select a joint point by pressing ‘J’ and clicking anywhere on the cage : the closest vertex to the click will be set as the joint point. The point will appear in blue.
- Select a bunch of points to move by clicking and dragging the mouse : click at approximately ¾ between the joint and the end of the limb you want to move in order to cover the whole limb. 
Dragging the mouse (while still pressing “S”) will rotate and stretch those selected points in the plane parallel to the image plane.
- Click anywhere on the mesh while pressing ‘J’ to unselect the joint point, and restart the procedure.
