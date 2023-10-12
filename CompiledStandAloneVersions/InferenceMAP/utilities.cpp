///*
// *
// *  InferenceMAP v1.0
// *  22/05/2014
// *
// *  Author: Mohamed El Beheiry, Physico-Chimie Curie, Institut Curie
// *  		Jean-Baptiste Masson, Physics of Biological Systems, Institut Pasteur
// *  Contact e-mail: mohamed.elbeheiry@gmail.com
// *  Copyright (c) 2014, Mohamed El Beheiry, Jean-Baptiste Masson, Institut Curie, Institut Pasteur
// *  All rights Reserved.
// *
// *  InferenceMAP is released under an "academic use only" licence.
// *  Details of which are provided in the xxx.doc file.
// *  Usage of InferenceMAP requires acceptance of this license.
// *
// *  User instructions for using InferenceMAP are provided in the InferenceMAP User Manual.
// *
// */
//
//#include "utilities.h"
//
//// siz             size of the N-dimensional matrix
//// N               the dimensions of the matrix
//// idx             index in linear format
//// sub             the output - subscript values written into an N-dimensional integer array
//// Example for a 2-D array of size [3, 5] : siz[] = {3, 5}; N = 2
//
//void ind2sub(int *siz, int N, int idx, int *sub) {
//	if (0)
//	{
//		int prod[4];
//		prod[0] = siz[3]*siz[2]*siz[1];
//		prod[1] = siz[3]*siz[2];
//		prod[2] = siz[3];
//		prod[3] = 1;
//		sub[0] = (int)floor(	(float)idx / prod[0]							);
//		sub[1] = (int)floor(	(float)(	idx % prod[0]	)/prod[1]				);
//		sub[2] = (int)floor( (float)( (	idx % prod[0]	)%prod[1]	)  / prod[2]);
//		sub[3] = 		( (	idx % prod[0]	)%prod[1]	)  % prod[2] ;
//	}
//	else
//	{
//		int *prod = new int [N];
//		for (int i = 0; i < N; i++)
//		{
//			prod[i] = 1;
//			for (int j = N-1; j > i; j--)
//				prod[i] *= siz[j];
//		}
//		//sub[0] = (int)floor( (float)idx / prod[0] );
//
//		for (int i = 0; i < N; i++)
//		{
//			sub[i] = idx ;
//			for (int j = 0; j < i ; j++)
//				sub[i] = sub[i] % prod[j];
//			sub[i] = (int)floor( (float)sub[i] / prod[i] );
//		}
//		//sub[N-1] = idx ;
//		//for (int j = 0; j < N -1; j++)
//		//	sub[N-1] = sub[N-1] % prod[j];
//
//		delete [] prod;
//	}
//}
//
//// siz             size of the N-dimensional matrix
//// N               the dimensions of the matrix
//// sub             subscripts stored in an N-dimensional array
//// return         (integer) the linear index corresponding to the subscript
//// Example for a 2-D array of size [3, 5] : siz[] = {3, 5}; N = 2
//
//int sub2ind(int *siz, int N, int *sub) {
//	int idx =	0;
//	if (0)
//		idx =		sub[0]*siz[3]*siz[2]*siz[1]
//				+	sub[1]*siz[3]*siz[2]
//				+	sub[2]*siz[3]
//				+	sub[3];
//	else
//	{
//		for (int i = 0; i < N; i++)
//		{
//			int prod = 1;
//			for (int j = N-1; j > i; j--)
//				prod *= siz[j];
//			idx += sub[i] * prod;
//		}
//	}
//	return idx;
//}
//
/////
////  Create an OpenCL context on the first available platform using
////  either a GPU or CPU depending on what is available.
////
//cl_context CreateContext()
//{
//    cl_int errNum;
//    cl_uint numPlatforms;
//    cl_platform_id firstPlatformId;
//    cl_context context = NULL;
//
//    // First, select an OpenCL platform to run on.  For this example, we
//    // simply choose the first available platform.  Normally, you would
//    // query for all available platforms and select the most appropriate one.
//    errNum = clGetPlatformIDs(1, &firstPlatformId, &numPlatforms);
//    if (errNum != CL_SUCCESS || numPlatforms <= 0)
//    {
//        fprintf(stderr,"Failed to find any OpenCL platforms.\n");
//        return NULL;
//    }
//
//    // Next, create an OpenCL context on the platform.  Attempt to
//    // create a GPU-based context, and if that fails, try to create
//    // a CPU-based context.
//    cl_context_properties contextProperties[] =
//    {
//        CL_CONTEXT_PLATFORM,
//        (cl_context_properties)firstPlatformId,
//        0
//    };
//    context = clCreateContextFromType(contextProperties, CL_DEVICE_TYPE_GPU,
//                                      NULL, NULL, &errNum);
//    if (errNum != CL_SUCCESS)
//    {
//        fprintf(stderr,"Could not create GPU context, trying CPU...\n");
//        context = clCreateContextFromType(contextProperties, CL_DEVICE_TYPE_CPU,
//                                          NULL, NULL, &errNum);
//        if (errNum != CL_SUCCESS)
//        {
//            fprintf(stderr,"Failed to create an OpenCL GPU or CPU context.\n");
//            return NULL;
//        }
//    }
//
//    return context;
//}
//
/////
////  Create a command queue on the first device available on the
////  context
////
//cl_command_queue CreateCommandQueue(cl_context context, cl_device_id *device)
//{
//    cl_int errNum;
//    cl_device_id *devices;
//    cl_command_queue commandQueue = NULL;
//    size_t deviceBufferSize = -1;
//
//    // First get the size of the devices buffer
//    errNum = clGetContextInfo(context, CL_CONTEXT_DEVICES, 0, NULL, &deviceBufferSize);
//    if (errNum != CL_SUCCESS)
//    {
//        fprintf(stderr,"Failed call to clGetContextInfo(...,GL_CONTEXT_DEVICES,...)\n");
//        return NULL;
//    }
//
//    if (deviceBufferSize <= 0)
//    {
//        fprintf(stderr,"No devices available.\n");
//        return NULL;
//    }
//
//    // Allocate memory for the devices buffer
//    devices = new cl_device_id[deviceBufferSize / sizeof(cl_device_id)];
//    errNum = clGetContextInfo(context, CL_CONTEXT_DEVICES, deviceBufferSize, devices, NULL);
//    if (errNum != CL_SUCCESS)
//    {
//        delete [] devices;
//        fprintf(stderr,"Failed to get device IDs\n");
//        return NULL;
//    }
//
//    // In this example, we just choose the first available device.  In a
//    // real program, you would likely use all available devices or choose
//    // the highest performance device based on OpenCL device queries
//    commandQueue = clCreateCommandQueue(context, devices[0], 0, NULL);
//    if (commandQueue == NULL)
//    {
//        delete [] devices;
//        fprintf(stderr,"Failed to create commandQueue for device 0\n");
//        return NULL;
//    }
//
//    *device = devices[0];
//    delete [] devices;
//    return commandQueue;
//}
//
/////
////  Create an OpenCL program from the kernel source file
////
//cl_program CreateProgram(cl_context context, cl_device_id device, const char* fileName)
//{
//    cl_int errNum;
//    cl_program program;
//
//    std::ifstream kernelFile(fileName, std::ios::in);
//    if (!kernelFile.is_open())
//    {
//    	fprintf(stderr,"Failed to open file for reading: %s\n",fileName);
//        return NULL;
//    }
//
//    std::ostringstream oss;
//    oss << kernelFile.rdbuf();
//
//    std::string srcStdStr = oss.str();
//    const char *srcStr = srcStdStr.c_str();
//    program = clCreateProgramWithSource(context, 1,
//                                        (const char**)&srcStr,
//                                        NULL, NULL);
//    if (program == NULL)
//    {
//    	fprintf(stderr,"Failed to create CL program from source.\n");
//        return NULL;
//    }
//
//    errNum = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
//    if (errNum != CL_SUCCESS)
//    {
//        // Determine the reason for the error
//        char buildLog[16384];
//        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG,
//                              sizeof(buildLog), buildLog, NULL);
//
//        fprintf(stderr,"Error in kernel: \n%s\n",buildLog);
//        clReleaseProgram(program);
//        return NULL;
//    }
//
//    return program;
//}
//
/////
////  Create memory objects used as the arguments to the kernel
////  The kernel takes three arguments: result (output), a (input),
////  and b (input)
////
//bool CreateMemObjects(cl_context context, cl_mem memObjects[3], float *a, float *b)
//{
////    memObjects[0] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
////                                   sizeof(float) * ARRAY_SIZE, a, NULL);
////    memObjects[1] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
////                                   sizeof(float) * ARRAY_SIZE, b, NULL);
////    memObjects[2] = clCreateBuffer(context, CL_MEM_READ_WRITE,
////                                   sizeof(float) * ARRAY_SIZE, NULL, NULL);
////
////    if (memObjects[0] == NULL || memObjects[1] == NULL || memObjects[2] == NULL)
////    {
////    	fprintf(stderr,"Error creating memory objects.\n");
////        return false;
////    }
//
//    return true;
//}
//
/////
////  Cleanup any created OpenCL resources
////
//void Cleanup(cl_context context, cl_command_queue commandQueue, cl_program program, cl_kernel kernel, cl_mem memObjects[3])
//{
//    for (int i = 0; i < 3; i++)
//    {
//        if (memObjects[i] != 0)
//            clReleaseMemObject(memObjects[i]);
//    }
//    if (commandQueue != 0)
//        clReleaseCommandQueue(commandQueue);
//
//    if (kernel != 0)
//        clReleaseKernel(kernel);
//
//    if (program != 0)
//        clReleaseProgram(program);
//
//    if (context != 0)
//        clReleaseContext(context);
//
//}
//
///***** CAMERA CLASS *****/
//
//GLfloat GetVectorLength(Vector3D * v) {
//	return (GLfloat)(sqrt((v->x)*(v->x)+(v->y)*(v->y)+(v->z)*(v->z)));
//}
//
//Vector3D Normalize3dVector(Vector3D v) {
//	Vector3D res;
//	float l = GetVectorLength(&v);
////	if (l == 0.0f) return NULL;
//	res.x = v.x / l;
//	res.y = v.y / l;
//	res.z = v.z / l;
//	return res;
//}
//
//Vector3D operator+ (Vector3D v, Vector3D u) {
//	Vector3D res;
//	res.x = v.x+u.x;
//	res.y = v.y+u.y;
//	res.z = v.z+u.z;
//	return res;
//}
//Vector3D operator- (Vector3D v, Vector3D u) {
//	Vector3D res;
//	res.x = v.x-u.x;
//	res.y = v.y-u.y;
//	res.z = v.z-u.z;
//	return res;
//}
//
//Vector3D operator* (Vector3D v, float r) {
//	Vector3D res;
//	res.x = v.x*r;
//	res.y = v.y*r;
//	res.z = v.z*r;
//	return res;
//}
//
//Vector3D CrossProduct (Vector3D * u, Vector3D * v) {
//	Vector3D resVector;
//	resVector.x = u->y*v->z - u->z*v->y;
//	resVector.y = u->z*v->x - u->x*v->z;
//	resVector.z = u->x*v->y - u->y*v->x;
//
//	return resVector;
//}
//
//float DotProduct (Vector3D * u, Vector3D * v) {
//	return u->x*v->x + u->y*v->y + u->z*v->z;
//}
//
//float operator* (Vector3D v, Vector3D u) {
//	return v.x*u.x+v.y*u.y+v.z*u.z;
//}
//
//Camera::Camera() {
//	//Init with standard OGL values:
//	Position = Vector3D (0.0, 0.0,	0.0);
//	ViewDir = Vector3D( 0.0, 0.0, -1.0);
//	RightVector = Vector3D (1.0, 0.0, 0.0);
//	UpVector = Vector3D (0.0, 1.0, 0.0);
//
//	//Only to be sure:
//	RotatedX = RotatedY = RotatedZ = 0.0;
//}
//
//void Camera::Move (Vector3D Direction) {
//	Position = Position + Direction;
//}
//
//void Camera::RotateX (GLfloat Angle) {
//	RotatedX += Angle;
//	//Rotate viewdir around the right vector:
//	ViewDir = Normalize3dVector(ViewDir*cos(Angle*PIdiv180)
//								+ UpVector*sin(Angle*PIdiv180));
//	//now compute the new UpVector (by cross product)
//	UpVector = CrossProduct(&ViewDir, &RightVector)*-1;
//}
//
//void Camera::RotateY (GLfloat Angle) {
//	RotatedY += Angle;
//	//Rotate viewdir around the up vector:
//	ViewDir = Normalize3dVector(ViewDir*cos(Angle*PIdiv180)
//								- RightVector*sin(Angle*PIdiv180));
//	//now compute the new RightVector (by cross product)
//	RightVector = CrossProduct(&ViewDir, &UpVector);
//}
//
//void Camera::RotateZ (GLfloat Angle) {
//	RotatedZ += Angle;
//	//Rotate viewdir around the right vector:
//	RightVector = Normalize3dVector(RightVector*cos(Angle*PIdiv180)
//								+ UpVector*sin(Angle*PIdiv180));
//	//now compute the new UpVector (by cross product)
//	UpVector = CrossProduct(&ViewDir, &RightVector)*-1;
//}
//
//void Camera::Render( void ) {
//
//	//The point at which the camera looks:
//	Vector3D ViewPoint = Position+ViewDir;
//
//	//as we know the up vector, we can easily use gluLookAt:
//	gluLookAt(	Position.x,Position.y,Position.z,
//				ViewPoint.x,ViewPoint.y,ViewPoint.z,
//				UpVector.x,UpVector.y,UpVector.z);
//
//}
//
//void Camera::MoveForward( GLfloat Distance ) {
//	Position = Position + (ViewDir*-Distance);
//}
//
//void Camera::StrafeRight ( GLfloat Distance ) {
//	Position = Position + (RightVector*Distance);
//}
//
//void Camera::MoveUpward( GLfloat Distance ) {
//	Position = Position + (UpVector*Distance);
//}
//
//void Camera::Reset( void ) {
//	Position = Vector3D (0.0, 0.0,	0.0);
//	ViewDir = Vector3D( 0.0, 0.0, -1.0);
//	RightVector = Vector3D (1.0, 0.0, 0.0);
//	UpVector = Vector3D (0.0, 1.0, 0.0);
//	RotatedX = 0.0;
//	RotatedY = 0.0;
//	RotatedZ = 0.0;
//}
