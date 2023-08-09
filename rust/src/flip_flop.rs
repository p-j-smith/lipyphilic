use numpy::ndarray::ArrayView1;

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
pub(crate) fn molecule_flip_flop(
    leaflets: ArrayView1<i8>,
    frame_cutoff: usize,
) -> (Vec<usize>, Vec<usize>, Vec<isize>, Vec<String>) {
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

    return (start_frames, end_frames, end_leaflets, success);
}
