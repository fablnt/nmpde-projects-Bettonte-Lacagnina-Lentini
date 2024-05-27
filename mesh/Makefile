# Define the Gmsh executable
GMSH = gmsh

# Define the target .msh file
TARGET = ellipse.msh

# Define the source .geo file
SOURCE = new_ellipse.geo

# Define the mesh size scaling factor (adjust this as needed)
MESH_SIZE_SCALE = 0.06

# Default rule to generate .msh from .geo
all: $(TARGET)

$(TARGET): $(SOURCE)
	$(GMSH) -2 $(SOURCE) -clscale $(MESH_SIZE_SCALE) -o $(TARGET)

# Clean rule to remove generated .msh file
clean:
	rm -f $(TARGET)
