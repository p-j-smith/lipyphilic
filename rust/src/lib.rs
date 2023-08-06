use numpy::borrow::PyReadonlyArray1;
use pyo3::types::{PyList, PyModule};
use pyo3::{pyfunction, pymodule, wrap_pyfunction, Py, PyResult, Python};

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

/// Find flip-flop events for a single molecule based on its leaflet membership
/// Leflet membership must be defined by i8 integers as follows:
/// -1: lower leaflet
/// 0: midplane
/// 1: upper leaflet
#[pyfunction]
fn molecule_flip_flop(
    leaflets: PyReadonlyArray1<i8>,
    frame_cutoff: i32,
) -> PyResult<(Py<PyList>, Py<PyList>, Py<PyList>, Py<PyList>)> {
    let start_frames: Vec<i32> = Vec::new();
    let end_frames: Vec<i32> = Vec::new();
    let end_leaflets: Vec<i32> = Vec::new();
    let success: Vec<String> = Vec::new();

    Python::with_gil(|py| {
        let py_start_frames: &PyList = PyList::new(py, &start_frames);
        let py_end_frames: &PyList = PyList::new(py, &end_frames);
        let py_end_leaflets: &PyList = PyList::new(py, &end_leaflets);
        let py_success: &PyList = PyList::new(py, &success);

        Ok((
            py_start_frames.into(),
            py_end_frames.into(),
            py_end_leaflets.into(),
            py_success.into(),
        ))
    })
}

/// A Python module implemented in Rust.
#[pymodule]
fn _lipyferrous(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    m.add_function(wrap_pyfunction!(molecule_flip_flop, m)?)?;
    Ok(())
}
