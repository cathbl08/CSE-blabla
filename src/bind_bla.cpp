#include <sstream>
#include <pybind11/pybind11.h>
#include <cstring>

#include "vector.hpp"
#include "matrix.hpp"
#include "matexpr.hpp"
#include "lapack_interface.hpp"

using namespace ASC_bla;
namespace py = pybind11;




PYBIND11_MODULE(bla, m) {
    m.doc() = "Basic linear algebra module"; // optional module docstring
    
    py::class_<Vector<double>> (m, "Vector")
      .def(py::init<size_t>(),
           py::arg("size"), "create vector of given size")
      .def("__len__", &Vector<double>::size,
           "return size of vector")
      
      .def("__setitem__", [](Vector<double> & self, int i, double v) {
        if (i < 0) i += self.size();
        if (i < 0 || i >= self.size()) throw py::index_error("vector index out of range");
        self(i) = v;
      })
      .def("__getitem__", [](Vector<double> & self, int i) { return self(i); })
      
      .def("__setitem__", [](Vector<double> & self, py::slice inds, double val)
      {
        size_t start, stop, step, n;
        if (!inds.compute(self.size(), &start, &stop, &step, &n))
          throw py::error_already_set();
        self.range(start, stop).slice(0,step) = val;
      })
      
      .def("__add__", [](Vector<double> & self, Vector<double> & other)
      { return Vector<double> (self+other); })

      .def("__rmul__", [](Vector<double> & self, double scal)
      { return Vector<double> (scal*self); })

      .def("__mul__", [](Vector<double> & self, double scal)
      { return Vector<double> (scal*self); })
      
      .def("__str__", [](const Vector<double> & self)
      {
        std::stringstream str;
        str << self;
        return str.str();
      })

     .def(py::pickle(
        [](Vector<double> & self) { // __getstate__
            /* return a tuple that fully encodes the state of the object */
          return py::make_tuple(self.size(),
                                py::bytes((char*)(void*)&self(0), self.size()*sizeof(double)));
        },
        [](py::tuple t) { // __setstate__
          if (t.size() != 2)
            throw std::runtime_error("should be a 2-tuple!");

          Vector<double> v(t[0].cast<size_t>());
          py::bytes mem = t[1].cast<py::bytes>();
          std::memcpy(&v(0), PYBIND11_BYTES_AS_STRING(mem.ptr()), v.size()*sizeof(double));
          return v;
        }));


      //MATRIX 

    py::class_<Matrix<double>> (m, "Matrix")
    .def(py::init<size_t, size_t>(),
         py::arg("rows"), py::arg("cols"), "create matrix of given size")

    .def_property_readonly("shape",
      [](const Matrix<double, RowMajor>& self) {
           return std::tuple(self.rows(), self.cols());
      })

    .def("__getitem__",
      [](Matrix<double, RowMajor> self, std::tuple<int, int> ind) {
           return self(std::get<0>(ind), std::get<1>(ind));
      })

    .def("__setitem__", [](Matrix<double> & self, py::tuple i, double m) {
        if (i.size() != 2)
          throw py::index_error("matrix index must be a 2-tuple");
        size_t row = i[0].cast<size_t>();
        size_t col = i[1].cast<size_t>();
        self(row, col) = m;
    })

    .def("__add__", [](Matrix<double> & self, Matrix<double> & other)
      { return Matrix<double> (self+other); })

    // scalar * Matrix and Matrix * scalar still fine:
    .def("__rmul__", [](Matrix<double> & self, double scal)
      { return Matrix<double, RowMajor> (scal * self); })

    .def("__mul__", [](Matrix<double> & self, double scal)
      { return Matrix<double, RowMajor> (self * scal); })

    // ------- THIS is the important change: Matrix * Matrix via ET -------
    .def("__mul__", [](const Matrix<double>& A, const Matrix<double>& B) {
        if (A.cols() != B.rows())
            throw std::runtime_error("Incompatible shapes");

        Matrix<double> C(A.rows(), B.cols());

        // Views (these inherit from MatExpr, so Av*Bv builds MultMatExpr)
        MatrixView<double> Cv = C;
        MatrixView<double> Av = const_cast<Matrix<double>&>(A);
        MatrixView<double> Bv = const_cast<Matrix<double>&>(B);

        Cv = Av * Bv;   // <- expression templates used here
        return C;
    })

    .def("__matmul__", [](const Matrix<double>& A, const Matrix<double>& B) {
        // you can choose: either call A*B again, or use LAPACK here, etc.
        return A * B;
    })

    .def("__str__", [](const Matrix<double> & self)
      {
        std::stringstream str;
        str << self;
        return str.str();
      })

    .def(py::pickle(
        [](Matrix<double> & self) {
          return py::make_tuple(self.rows(), self.cols(),
                                py::bytes((char*)(void*)&self(0, 0),
                                          self.rows()*self.cols()*sizeof(double)));
        },
        [](py::tuple t) {
          if (t.size() != 2)
            throw std::runtime_error("should be a 2-tuple!");

          Matrix<double> m(t[0].cast<size_t>());
          py::bytes mem = t[1].cast<py::bytes>();
          std::memcpy(&m(0, 0), PYBIND11_BYTES_AS_STRING(mem.ptr()),
                      m.rows()*m.cols()*sizeof(double));
          return m;
        }))
    ;

    m.def("matmul_lapack", [](const Matrix<double, RowMajor>& A,
                          const Matrix<double, RowMajor>& B)
    {
        if (A.cols() != B.rows())
            throw std::runtime_error("Incompatible shapes");

        Matrix<double, RowMajor> C(A.rows(), B.cols());

        MatrixView<double, RowMajor> Av = A;
        MatrixView<double, RowMajor> Bv = B;
        MatrixView<double, RowMajor> Cv = C;

        Cv = (Av * Bv) | Lapack;  // calls dgemm
        return C;
    }, "Matrix-matrix multiply using LAPACK (dgemm)");
}