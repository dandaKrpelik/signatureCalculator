
#include "graph.h"
#include "system.h"
#include "pChrono.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>

using namespace std;

PYBIND11_MODULE(pyIPSIG, m)
{
	m.doc()="For calculating (imprecise) system survival signature from reliability block diagrams.";
	
	namespace py = pybind11;
	using namespace pybind11::literals;
  
	py::class_<Graph>(m, "Graph")
		.def(py::init< vector<vector<size_t> > >())
		.def_readonly("N", &Graph::_N)
		.def_readonly("E", &Graph::E)
		;

	py::class_<Branch>(m, "Branch")
		.def(py::init< vector<bool>, posvec, posvec >())
		.def_readonly("mask", &Branch::mask)
		.def_readonly("ones", &Branch::ones)
		.def_readonly("minpath", &Branch::minpath)
		;
		

	py::class_<System>(m, "System")
		.def(py::init< Graph*, vector<size_t> >())
		.def_readonly("N", &System::N)
		.def_readonly("K", &System::K)
		.def_readonly("M", &System::M)
		.def_readonly("C", &System::C)
		.def_readonly("g", &System::g)
		.def_readonly("branch", &System::branch)
		.def_readonly("sig1", &System::sig1)
		.def_readonly("sig0", &System::sig0)
		.def_readonly("mul", &System::mul)
		.def_readonly("binCo", &System::binCo)
		.def("pos2index", &System::pos2index)
		.def("index2pos", &System::index2pos)
		.def("iterate", &System::iterate)
		.def("SIG", &System::SIGp)
		.def("calcSIG", &System::calcSIG)
		.def("eval", &System::eval)
		.def_readonly("LowSig", &System::LowSig)
		.def_readonly("HiSig", &System::HiSig)
		.def_readonly("validSig", &System::validSig)
		;
  
    
	m.def("helloWorld",&helloWorld,"a"_a=1)
		.def("chStart", &chStart)
		.def("chEnd", &chEnd)
		.def("minPath", &minPath)
		.def("minPath_plain", &minPath_plain)
		.def("binomialCoeff", &binomialCoeff)
	;
    
}
