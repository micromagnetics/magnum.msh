%module(package="magnummsh") cpp

%{
#include "Mesher.h"
%}

// Handle NumPy and Dolfin Arrays
%{
#include <numpy/arrayobject.h> 
#define SWIG_SHARED_PTR_QNAMESPACE boost
%}

%include <exception.i>

%init%{
import_array();
%}

%import "typemaps/numpy.i"
%import "typemaps/array.i"
%import "typemaps/std_map.i"
%import "typemaps/std_vector.i"
%import "typemaps/primitives.i"

// Handle Strings
%include <std_string.i>

// Handle Shared Pointers
%include <std_shared_ptr.i>
%shared_ptr(dolfin::FunctionSpace)
%shared_ptr(dolfin::GenericTensor)
%shared_ptr(dolfin::GenericVector)
%shared_ptr(dolfin::GenericMatrix)
%shared_ptr(dolfin::GenericFunction)
%shared_ptr(dolfin::Mesh)
%shared_ptr(dolfin::SubDomain)

// Include headers
%include "Mesher.h"

// vim:ft=cpp:
