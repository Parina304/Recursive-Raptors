#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPolyDataWriter.h>

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " input.vtk output.vtk" << std::endl;
        return EXIT_FAILURE;
    }

    // Read the VTK file
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(argv[1]);
    reader->Update();

    vtkSmartPointer<vtkPolyData> mesh = reader->GetOutput();

    // Create a transformation
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    
    // Translation (move the object)
    transform->Translate(10.0, 5.0, 0.0);  // Example: Translate by (10,5,0)

    // Rotation (rotate around an axis)
    transform->RotateX(30.0);  // Rotate 30 degrees around X-axis
    transform->RotateY(45.0);  // Rotate 45 degrees around Y-axis
    transform->RotateZ(60.0);  // Rotate 60 degrees around Z-axis

    // Scaling (resize the object)
    transform->Scale(1.5, 1.5, 1.5);  // Scale by 1.5x in all directions

    // Reflection (mirror across an axis)
    // Use -1 scaling in any axis to reflect
    // transform->Scale(-1.0, 1.0, 1.0);  // Reflect across the YZ-plane (mirror along X-axis)

    // Apply the transformation
    vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transformFilter->SetInputData(mesh);
    transformFilter->SetTransform(transform);
    transformFilter->Update();

    // Write the transformed mesh to a new VTK file
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(argv[2]);
    writer->SetInputData(transformFilter->GetOutput());
    writer->Write();

    std::cout << "Transformations applied and saved to " << argv[2] << std::endl;

    return EXIT_SUCCESS;
}
