use numpy::borrow::PyReadonlyArray1;
use numpy::ndarray::ArrayView1;
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
    frame_cutoff: usize,
) -> PyResult<(Py<PyList>, Py<PyList>, Py<PyList>, Py<PyList>)> {
    let leaflets: ArrayView1<i8> = leaflets.as_array();
    let mut start_frames: Vec<usize> = Vec::new();
    let mut end_frames: Vec<usize> = Vec::new();
    let mut end_leaflets: Vec<isize> = Vec::new();
    let mut success: Vec<String> = Vec::new();

    let mut current_leaflet: i8 = leaflets[0]; // leaflet we're currently in
    let mut apposing_leaflet: i8 = current_leaflet * -1;

    let mut event_start: Option<usize> = None;
    let mut event_stop: Option<usize> = None;
    let mut event_success: Option<bool> = None;

    let mut n_left: usize = 0; // number of consecutive frames lipid has left its current leaflet
    let mut n_apposing: usize = 0; // number of consecutive frames lipid as been in apposing leaflet
    let mut n_returned: usize = 0; // number of consecutive frames lipid has returned to current leaflet since event started

    let n_frames: usize = leaflets.shape()[0];

    for frame in 1..n_frames {
        let leaflet: i8 = leaflets[frame];

        if event_start.is_none() {
            // Check if an event has started and reset the counter if so
            n_left = if leaflet == current_leaflet {
                0
            } else {
                n_left + 1
            };
            event_start = if n_left == frame_cutoff {
                Some(frame - frame_cutoff)
            } else {
                None
            };
            n_left = if event_start.is_none() { n_left } else { 0 };
        }

        n_apposing = if leaflet != apposing_leaflet {
            0
        } else {
            n_apposing + 1
        };

        if event_start.is_none() {
            continue;
        }

        n_returned = if leaflet != current_leaflet {
            0
        } else {
            n_returned + 1
        };

        // Failed event
        if n_returned == frame_cutoff {
            event_stop = Some(frame - (frame_cutoff - 1));
            event_success = Some(false);
        }

        // Successful event
        if n_apposing == frame_cutoff {
            event_stop = Some(frame - (frame_cutoff - 1));
            event_success = Some(true);
        }

        // Event hasn't ended
        if event_stop.is_none() {
            continue;
        }

        // Event has ended! Swap leaflets if necessary
        if let Some(succeeded) = event_success {
            if succeeded {
                (current_leaflet, apposing_leaflet) = (apposing_leaflet, current_leaflet);
            }
        }

        start_frames.push(event_start.unwrap());
        end_frames.push(event_stop.unwrap());
        end_leaflets.push(current_leaflet as isize);
        success.push(if event_success.unwrap() {
            "Success".to_string()
        } else {
            "Fail".to_string()
        });

        // Reset all counters
        event_start = None;
        event_stop = None;
        event_success = None;

        n_left = 0;
        n_apposing = 0;
        n_returned = 0;
    }

    if let Some(start) = event_start {
        start_frames.push(start);
        end_frames.push(n_frames - 1);
        end_leaflets.push(leaflets[n_frames - 1] as isize);
        success.push("Ongoing".to_string());
    }

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
