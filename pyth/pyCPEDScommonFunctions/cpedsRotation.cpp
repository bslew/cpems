/*!
 * \file pycpeds-math.cpp
 *
 *  Created on: Aug 29, 2018
 *      Author: blew
 */



#include <Python.h>
#include <ndarraytypes.h>
#include <numpy/arrayobject.h>
#include "cpeds-math.h"
#include "cpeds-point3d.h"

static PyObject *rotX(PyObject *self, PyObject *args) {
	double x,y,z,ang;
	if (!PyArg_ParseTuple(args, "dddd", &x,&y,&z,&ang)) {
		return NULL;
	}
	//  printf("x: %lf, y: %lf, z: %lf, ang: %lf\n",x,y,z,ang);
	cpedsPoint3D p(x,y,z);
	p.Rx(ang);
	
	//  printf("x: %lf, y: %lf, z: %lf, ang: %lf\n",x,y,z,ang);
	//  PyObject* l= PyList_New(3);
	//  PyObject* px=PyFloat_FromDouble(p.x());
	//  PyObject* py=PyFloat_FromDouble(p.y());
	//  PyObject* pz=PyFloat_FromDouble(p.z());
	//  
	//  PyList_SetItem(l,0,px);
	//  PyList_SetItem(l,1,py);
	//  PyList_SetItem(l,2,pz);
	//  return l;
	
	return Py_BuildValue("ddd",p.x(),p.y(),p.z());
}

static PyObject *rotY(PyObject *self, PyObject *args) {
	double x,y,z,ang;
	if (!PyArg_ParseTuple(args, "dddd", &x,&y,&z,&ang)) {
		return NULL;
	}
	cpedsPoint3D p(x,y,z);
	p.Ry(ang);
	return Py_BuildValue("ddd",p.x(),p.y(),p.z());
}

static PyObject *rotZ(PyObject *self, PyObject *args) {
	double x,y,z,ang;
	if (!PyArg_ParseTuple(args, "dddd", &x,&y,&z,&ang)) {
		return NULL;
	}
	cpedsPoint3D p(x,y,z);
	p.Rz(ang);
	return Py_BuildValue("ddd",p.x(),p.y(),p.z());
}
/* ******************************************************************************************** */
/*!
	\brief rotate set of points by Rx
	\details 
	@param
	@return

	\date Aug 29, 2018, 7:30:53 PM
 */
static PyObject *rotXarray(PyObject *self, PyObject *args) {
	PyObject *a=NULL;
	double ang;
	
	if (!PyArg_ParseTuple(args, "dO", &ang,&a)) return NULL;
	//	printf("ndim: %i\n",(int)PyArray_NDIM(arg1));
	
//	if (PyArray_NDIM(arg1) != 2 || PyArray_TYPE(arg1) != NPY_DOUBLE) {
//		PyErr_Format(PyExc_ValueError,
//				"a array is %d-dimensional or not of type double",
//				PyArray_NDIM(arg1));
//		return NULL;
//	}
	if (PyArray_NDIM(a) != 2) {
		PyErr_Format(PyExc_ValueError,
				"array is %d-dimensional (should be 2-dimensional)",
				PyArray_NDIM(a));
		return NULL;
	}
	
	int Nx=0,Ny=0;
	double **aC;
	double *ap;
	Nx = PyArray_DIM(a,0); 
	Ny = PyArray_DIM(a,1);
//	printf("nx: %i, ny: %i\n",Nx,Ny);
	
	if (PyArray_DIM(a,1) <3) {
		PyErr_Format(PyExc_ValueError,
				"array is too small (should have at least 3 columns)");
		return NULL;
	}
	
	/* allocate help array for creating a double pointer: */
	aC = (double **) malloc(Nx*sizeof(double*));
	ap = (double *) PyArray_DATA(a);
	for (long i = 0; i < Nx; i++) { aC[i] = &(ap[i*Ny]); }
	
//	npy_intp app_dims[2];
//	int typenum = PyArray_Type(a);
//	PyArray_Descr *descr= PyArray_DescrFromType(NPY_DOUBLE);
//	exit(0);
//	PyArray_AsCArray(&a, (void*)&app, app_dims, 2, descr);
//	PyArray_AsCArray(&a, (void*)&app, app_dims, 2, NPY_DOUBLE);
	//	PyArray_Free(a, (void*) &app);
	//	if ((int)PyArray_NDIM(arg1)!=2) return NULL;
	
	/*
	PyObject *arr1=PyArray_FROM_OTF(arg1,NPY_DOUBLE,NPY_IN_ARRAY);
	exit(0);
	if (arr1 == NULL) {
		Py_XDECREF(arr1);
		return NULL;
	}
	 */
	
	/* How many data points are there? */
	//	int Nx = (int)PyArray_DIM(arr1, 0);
	//	int Ny = (int)PyArray_DIM(arr1, 1);
	//	printf("Nx: %i\n",Nx);
	//	printf("Ny: %i\n",Ny);
	
	/* Get pointers to the data as C-types. */
//		double *mat1    = (double*)PyArray_DATA(a);
//		printf("mat1[0]: %lE\n",mat1[0]);

	npy_intp a_dims[2]; a_dims[0] = Nx; a_dims[1] = Ny;
	PyArrayObject * aout = (PyArrayObject *) PyArray_SimpleNew(2, a_dims, NPY_DOUBLE);

	if (aout == NULL) {
		printf("creating %dx%d array failed\n",
				(int) a_dims[0], (int) a_dims[1]);
		return NULL; /* PyArray_SimpleNew raised an exception */
	}

	for (long i = 0; i < Nx; i++) {
//		for (long j = 0; j < Ny; j++) {
//			double* a_ij = (double *) PyArray_GETPTR2(a, i, j);
//			printf("a_ij: %lf\n",*a_ij);
//			printf("(%li, %li)=%lf\n",i,j,(double)aC[i][j]);
			cpedsPoint3D p(aC[i][0],aC[i][1],aC[i][2]);

			p.print_point("point");
			p.Rx(ang);
			p.print_point("rotated");
			
			double* x=(double*)PyArray_GETPTR2(aout, i, 0);
			double* y=(double*)PyArray_GETPTR2(aout, i, 1);
			double* z=(double*)PyArray_GETPTR2(aout, i, 2);
			*x=p.x();
			*y=p.y();
			*z=p.z();
//		}
	}
	
	free(aC);
	/* Clean up. */
	//	Py_DECREF(arr1);
	
	//	if (value < 0.0) {
	//		PyErr_SetString(PyExc_RuntimeError,
	//				"Chi-squared returned an impossible value.");
	//		return NULL;
	//	}
	
	/* Build the output tuple */
	//	PyObject *ret = Py_BuildValue("d", value);
	//	return ret;	
	
	
	//	return Py_BuildValue("ddd",p.x(),p.y(),p.z());
//	return NULL;
	return PyArray_Return(aout);
//	return Py_BuildValue("");
}




static PyMethodDef cpedsRotation_methods[] = {
		{"rotX", rotX, METH_VARARGS, "Returns Rx rotated x,y,z point coordinates"},
		{"rotY", rotY, METH_VARARGS, "Returns Ry rotated x,y,z point coordinates"},
		{"rotZ", rotZ, METH_VARARGS, "Returns Rz rotated x,y,z point coordinates"},
		{"rotXarray", rotXarray, METH_VARARGS, "rotates provided array of points by given angle"},
		{NULL, NULL, 0, NULL},
};

// python3 implementation
/*
static struct PyModuleDef pycpeds_math_definition = {
    PyModuleDef_HEAD_INIT,
    "example",
    "example module containing pants() function",
    -1,
    pycpeds_math_methods,
};

PyMODINIT_FUNC PyInit_example(void) {
  Py_Initialize();
  PyObject *m = PyModule_Create(&pycpeds_math_definition);

  return m;
}
 */

// python2 implementation
PyMODINIT_FUNC
initcpedsRotation(void)
{
	(void) Py_InitModule("cpedsRotation", cpedsRotation_methods);
	import_array();
}
