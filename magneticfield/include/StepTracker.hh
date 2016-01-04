// Nystrom stepper implemenation by Jason Suagee
//  Supervision / code review: John Apostolakis
//
// Sponsored by Google in Google Summer of Code 2015.
//
// First version: 27 May 2015
//
// This code is made available subject to the Geant4 license, a copy of
// which is available at
//   http://geant4.org/license

#ifndef MAGNETICFIELD_INCLUDE_STEPTRACKER_HH_
#define MAGNETICFIELD_INCLUDE_STEPTRACKER_HH_

#include "G4Types.hh"
#include "G4FieldTrack.hh"

#ifndef G4MAGIntegratorSTEPPER
#include "G4MagIntegratorStepper.hh"
#else
class G4MagIntegratorStepper;
#endif

#include "isTracking.hh"

#include <vector>
using namespace std;

#ifdef INTENDED_FOR_ERROR_BY_STEPPER_PROGRAM
#define BUFFER_COLUMN_LEN 28 // room for start point and end point of each step
// plus time/arclength entries for each.

#define ENDPOINT_BASE_INDEX 14
#define POSITION_SLOT 2
#define MOMENTUM_SLOT 5
#define RHS_SLOT 8
#define NUMBER_RHS_VARIABLES 6
#else
#define BUFFER_COLUMN_LEN 22  // room for start point and end point of each step
                              // plus time/arclength entries for each.
#define ENDPOINT_BASE_INDEX 11
#define POSITION_SLOT 2
#define MOMENTUM_SLOT 5
#define RHS_SLOT 8
#define NUMBER_RHS_VARIABLES 3
#endif

class StepTracker {
public:
    StepTracker();
    virtual ~StepTracker();

    // Initialization methods:
    void initialize_StepTracker(G4FieldTrack *initialFieldTrack);

    // Set the relativistic mass of the particle. This is only used
    // to calculate the starting velocity though. This does affect the
    // recorded times at each step.
    inline void set_mass(G4double mass_of_particle) { mass = mass_of_particle; }

    // Most important methods:
    void record_if_post_intersection_point( G4double passed_curve_length );

    virtual void RecordResultOfStepper(G4double yIn0[], G4double dydx0[],
                                       G4double yIn1[], G4double dydx1[],
                                       G4double arclength_to_add,
                                       G4int no_function_calls = -1);

    virtual void outputBuffer(char *outfile_name,
                              char *meta_outfile_name,
                              char *no_function_calls_outfile_name = 0,
                              char *no_function_calls_overshoot_filename = 0,
                              char *indices_intersection_pts_filename = 0,
                              char *overshoot_outfilename = 0);

    // Methods for communication and control handoff to StepTracker object
    // from within G4ChordFinder and G4MagIntegratorDriver objects:

    inline void set_within_AdvanceChordLimited(G4bool status) { within_AdvanceChordLimited = status; }
    void update_time_arclength(G4double time_to_add, G4double arclength_to_add);

    // alert StepTracker that the next step should be stored. This is to avoid storing
    // steps taking within DistChord() by an aux stepper, or extra steps involved with
    // using Richardson Extrapolation. (Currently Richardson Extrapolation is not an
    // issue because the two half steps are performed before the real step, and thus are
    // overwritten by the real step. However it's probably best to just be careful.)
    inline void ArmTracker() { armed = true; }

    inline void UnArmTracker() { armed = false; }
    inline void set_last_time_val_was_accepted(G4bool val) { last_time_val_was_accepted = val; }

    // Other methods

    inline G4bool get_within_AdvanceChordLimited() { return within_AdvanceChordLimited; }

    // This is used by ChordFinder and MagIntDriver to perform
    // the conversions from momentum coordinates to velocity
    // coordinates. It's needed currently because neither
    // ChordFinder or MagIntDriver do any record keeping.
    G4double last_velocity();

    inline bool isArmed() { return armed; }
    inline G4int getBufferLength() { return buffer.size(); }
    inline void add_to_num_other_function_calls_not_to_count(G4int k) { no_other_function_calls_not_to_count += k; }

private:
    G4int no_other_function_calls_not_to_count;

    G4bool within_AdvanceChordLimited;
    G4bool last_time_val_was_accepted, armed;

    G4double last_curve_length, last_time_length;
    G4double mass; // This is relativistic mass
    G4double first_velocity; // For when we need a velocity, but there are no items in the buffer.

    // Buffers used for output:
    vector<vector<G4double> > buffer;
    vector<G4int> no_function_calls_buffer;
    vector<G4int> indices_of_intersection_points;
    vector<vector<G4double> > overshoot_buffer;
    vector<G4int> no_function_calls_overshoot_buffer;
};

#endif /* MAGNETICFIELD_INCLUDE_STEPTRACKER_HH_ */
