use numpy::borrow::PyReadonlyArray1;
use numpy::ndarray::ArrayView1;
use pyo3::types::{PyList, PyModule};
use pyo3::{pyfunction, pymodule, wrap_pyfunction, Py, PyResult, Python};

pub(crate) mod flip_flop;

/// Find flip-flop events for a single molecule based on its leaflet membership
/// Leflet membership must be defined by i8 integers as follows:
/// -1: lower leaflet
/// 0: midplane
/// 1: upper leaflet
///
/// # Parameters
///
/// * `leaflets` - array of leaflet membership over time
/// * `frame_cutoff` - number of consecutive frames a molecule must
/// enter or leave a leaflet for an event to start or end
#[pyfunction]
fn molecule_flip_flop(
    leaflets: PyReadonlyArray1<i8>,
    frame_cutoff: usize,
) -> PyResult<(Py<PyList>, Py<PyList>, Py<PyList>, Py<PyList>)> {
    let leaflets: ArrayView1<i8> = leaflets.as_array();

    let (start_frames, end_frames, end_leaflets, success) =
        flip_flop::molecule_flip_flop(leaflets, frame_cutoff);

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
fn _lipyferric(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("version", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(molecule_flip_flop, m)?)?;
    Ok(())
}
