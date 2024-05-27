# Define the Gmsh executable
GMSH = gmsh

# Define the target .msh files
TARGET_1D = mesh1D.msh
TARGET_2D = ellipse.msh
TARGET_3D = brain-h3.0.msh

# Define the source .geo files
SOURCE_2D = new_ellipse.geo
SOURCE_1D = mesh1D.geo
SOURCE_3D = brain-h3.0.geo

# Define the mesh size scaling factors
MESH_SIZE_SCALE_1D = 0.01
MESH_SIZE_SCALE_2D = 0.06
MESH_SIZE_SCALE_3D = 0.1

# Default rule to generate .msh from .geo
all: $(TARGET_2D) $(TARGET_1D) $(TARGET_3D)

$(TARGET_2D): $(SOURCE_2D)
	$(GMSH) -2 $(SOURCE_2D) -clscale $(MESH_SIZE_SCALE_2D) -o $(TARGET_2D)

$(TARGET_1D): $(SOURCE_1D)
	$(GMSH) -1 $(SOURCE_1D) -clscale $(MESH_SIZE_SCALE_1D) -o $(TARGET_1D)

$(TARGET_3D): $(SOURCE_3D)
	$(GMSH) -3 $(SOURCE_3D) -clscale $(MESH_SIZE_SCALE_3D) -o $(TARGET_3D)

# Clean rule to remove generated .msh files
clean:
	rm -f $(TARGET_2D) $(TARGET_1D) $(TARGET_3D)
