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

#include "StepTracker.hh"
#include <fstream>
#include <vector>
using namespace std;
#include <assert.h>

#include "G4ThreeVector.hh"

#ifndef G4CACHED_MAGNETIC_FIELD_DEF
#include "G4CachedMagneticField.hh"
#endif

// For temporary testing with any future Nystrom stepper that rely upon an auxillary
// stepper for DistChord(). Use to profile how many function evaluations
// would be used by the stepper if interpolation methods were being used
// inside of DistChord() instead of an auxillary stepper.
#define IGNORING_DIST_CHORD_FUNCTION_CALLS


#define TOL 0.00000001


//StepTracker::StepTracker(G4double beginning[BUFFER_COLUMN_LEN]) {
   // beginning must be of the form
   // {time, arclength, position[1..3], momentum[1..3], RHS[1..3]}

   // Set this later in set_first_velocity() called from G4ChordFinder::AccurateAdvance()
   //first_velocity = G4ThreeVector( beginning[MOMENTUM_SLOT], beginning[MOMENTUM_SLOT + 1], beginning[MOMENTUM_SLOT + 2] ).mag();


StepTracker::StepTracker() {

   last_time_val_was_accepted = true; // There was no last time val at time of creation.
   armed = false;

   thrown_away_steps = 0;

   no_function_calls_buffer.push_back( 0 );
   no_function_calls_used_by_DistChord = 0;

   last_curve_length = last_time_length = 0.;

   within_AdvanceChordLimited = false;
}

StepTracker::~StepTracker() {
}


void StepTracker::initialize_StepTracker(G4FieldTrack *initialFieldTrack) {

   if (getBufferLength() == 0) {
      //G4ThreeVector Position = initialFieldTrack -> GetPosition();
      G4ThreeVector Momentum = initialFieldTrack->GetMomentum();
      first_velocity = Momentum.mag() / mass; // mass here is mass at rest.
   }
}

void StepTracker::record_if_post_intersection_point( G4FieldTrack& possible_post_intersection_point, //TODO: take this argument out, it isn't used
                                                     G4double passed_curve_length ) {

   assert( last_time_val_was_accepted );

   //cout << passed_curve_length << ", " << last_curve_length << endl;

   if ( passed_curve_length < last_curve_length ) {

      // passed_curve_length is the curve_length of the prev
      // Here is where we do something to record in the buffer that we had an intersection pt.

      // Note: this does not treat the case where the last point returned by AdvanceChordLimited()
      // was actually an intersection point (exactly). We might want to also consider the used
      // number of function evaluations to distinguish if there was no intersection in
      // PropagatorInField::ComputeStep(). If there was the code dealing with finding the intersection
      // points uses function evaluations (or at least function calls) (Is this true??)


      G4double last_length;
      // have to loop and throw away points until we have found pre-intersection point:

      vector<G4double> buffer_vector;

      do {
         //assert(buffer.size() > 0);

         last_length = buffer.back().at(ENDPOINT_BASE_INDEX + 1);

         buffer_vector = buffer.back();
         //temp_overshoot_buffer.push_back( buffer_vector ); // temp_overshoot_buffer will be in reverse order.

         buffer.pop_back(); // Get rid of overshoot point (overshot the intersection point).
         no_function_calls_buffer.pop_back(); // Similar as line above.
      } while ( buffer.size() > 0 &&
            passed_curve_length <= buffer.back().at(ENDPOINT_BASE_INDEX + 1) ); // TODO: check if this should be a < instead of a <=

      overshoot_buffer.push_back( buffer_vector );
      no_function_calls_overshoot_buffer.push_back( no_function_calls_buffer.back() );
      //no_function_calls_buffer.pop_back();


      G4double velocity;
      //differences_of_intersection_points.push_back( last_length - passed_curve_length );
      if (buffer.size() > 0) {
      last_time_length = buffer.back().at(ENDPOINT_BASE_INDEX + 0)                                 // Take the last pre-intesection time value
                           + ( passed_curve_length - buffer.back().at(ENDPOINT_BASE_INDEX + 1) )   // then add the remaining time (scaled arc length)
                              / last_velocity() ;
      }
      else {
         // This is what would have been the velocity, but we threw it away in the do-while loop
         velocity = G4ThreeVector( overshoot_buffer.back().at(MOMENTUM_SLOT + 0),
                                   overshoot_buffer.back().at(MOMENTUM_SLOT + 1),
                                   overshoot_buffer.back().at(MOMENTUM_SLOT + 2) ).mag();

         last_time_length = passed_curve_length / velocity; // We must have been at the beginning of the integration run.
      }

      last_curve_length = passed_curve_length; // Very important!

      if ( indices_of_intersection_points.size() > 0 ){
         if (indices_of_intersection_points.back() != getBufferLength() ) {
            assert( indices_of_intersection_points.back() < getBufferLength() );
            indices_of_intersection_points.push_back( getBufferLength() );
         }
      }
      else {
         //cout << passed_curve_length << ", " << passed_curve_length - buffer.back().at(ENDPOINT_BASE_INDEX + 1) << endl;
         indices_of_intersection_points.push_back( getBufferLength() );

      }
   }
   else {

      if (! (passed_curve_length < last_curve_length + TOL) ) {

         assert( passed_curve_length == last_curve_length );
      }

      //cout << "passed_curve_length: " << passed_curve_length << ",  last_curve_length: " << last_curve_length << endl;
   }
   // Either way we don't want to erase the last row of buffer data.
   last_time_val_was_accepted = true;
}

void StepTracker::outputBuffer(char *outfile_name,
                               char *meta_outfile_name,
                               char *no_function_calls_outfile_name,
                               char *no_function_calls_overshoot_filename,
                               char *indices_intersection_pts_filename,
                               //char *differences_of_intersection_points_filename
                               char * overshoot_outfilename) {

   ofstream meta_outfile(meta_outfile_name, ios::out);
   meta_outfile << getBufferLength() << endl << overshoot_buffer.size() << endl;
   meta_outfile.close();

   ofstream outfile;
   outfile.open(outfile_name, ios::binary | ios::out);

   for (int i = 0; i < getBufferLength(); i ++) {
      for (int j = 0; j < BUFFER_COLUMN_LEN; j ++){
         outfile.write( reinterpret_cast<char*>( &(buffer[i][j]) ), sizeof( G4double ) );
      }
   }

   outfile.close();

   if (no_function_calls_outfile_name != 0) {
      ofstream no_function_calls_outfile( no_function_calls_outfile_name, ios::out );
      for (int i = 0; i < getBufferLength(); i ++)
         no_function_calls_outfile << no_function_calls_buffer.at(i) << endl;

      no_function_calls_outfile.close();
   }

   if (indices_intersection_pts_filename != 0) {
      ofstream indices_intersection_pts_outfile( indices_intersection_pts_filename, ios::out );
      for (uint i = 0; i < indices_of_intersection_points.size(); i ++)
         indices_intersection_pts_outfile << indices_of_intersection_points[i] << endl;

      indices_intersection_pts_outfile.close();
   }

   if ( overshoot_outfilename != 0 ) {

      ofstream overshoot_outfile;
      overshoot_outfile.open(overshoot_outfilename, ios::binary | ios::out);

      for (int i = 0; i < overshoot_buffer.size(); i ++) {
         for (int j = 0; j < BUFFER_COLUMN_LEN; j ++){
            overshoot_outfile.write( reinterpret_cast<char*>( &(overshoot_buffer[i][j]) ),
                                                               sizeof( G4double ) );

         }
      }
      overshoot_outfile.close();
   }

   if (no_function_calls_overshoot_filename != 0) {
      ofstream no_function_calls_overshoot_file;
      no_function_calls_overshoot_file.open(no_function_calls_overshoot_filename, ios::binary | ios::out);

      for (int i = 0; i < no_function_calls_overshoot_buffer.size(); i ++) {
         no_function_calls_overshoot_file.write(
               reinterpret_cast<char*>( &(no_function_calls_overshoot_buffer[i]) ), sizeof( G4int ) );
      }
      no_function_calls_overshoot_file.close();
   }
}

void StepTracker::RecordResultOfStepper( G4double yIn0[], G4double dydx0[],
                                         G4double yIn1[], G4double dydx1[],
                                         G4double step,
                                         G4int no_function_calls) {
   // Time is stored in first component.

   if ( ! within_AdvanceChordLimited ) {

#ifdef DEBUG_TRACKING
      cout << "StepTracker::RecordResultOfStepper().  Not within AdvancedChordLimited()" << endl;
#endif
      return;
   }

   G4int last_index;

   // Create new space for next time/position/momentum values:
   if ( last_time_val_was_accepted ) {

      buffer.push_back( vector<G4double> (BUFFER_COLUMN_LEN) );

      last_index = buffer.size() - 1;


      // Copy last endpoint time/arclength values into new Initial point
      // time/arclength slots:

      buffer[last_index][0] = last_time_length;
      buffer[last_index][1] = last_curve_length;

      if (no_function_calls != -1)  // Never happens though

#ifdef IGNORING_DIST_CHORD_FUNCTION_CALLS

         no_function_calls_buffer.push_back( no_function_calls - no_function_calls_used_by_DistChord);
#else
         no_function_calls_buffer.push_back( no_function_calls );
#endif

      last_time_val_was_accepted = false;
   }
   // Otherwise, we will just write over the last position/momentum values
   else {
      thrown_away_steps ++;
   }


   last_index = buffer.size() - 1;

   // Copy over initial point values:
   for (int i = 0; i < 6; i ++) {
      buffer[last_index][i + POSITION_SLOT] = yIn0[i];
   }
   for (int i = 0; i < NUMBER_RHS_VARIABLES; i ++) {
      buffer[last_index][i + RHS_SLOT] = dydx0[i + 3];
   }
   // Now copy over endpoint values:
   for (int i = 0; i < 6; i ++) {
      buffer[last_index][ENDPOINT_BASE_INDEX + POSITION_SLOT + i] = yIn1[i];
   }
   for (int i = 0; i < NUMBER_RHS_VARIABLES; i ++) {
      buffer[last_index][ENDPOINT_BASE_INDEX + RHS_SLOT + i] = dydx1[i + 3];
   }

   buffer[last_index][ENDPOINT_BASE_INDEX + 0] = last_time_length + step / last_velocity();
   buffer[last_index][ENDPOINT_BASE_INDEX + 1] = last_curve_length + step;
}


G4double StepTracker::last_velocity() { // Might want to change this to use the endpoints velocity values?
   G4int last_index = getBufferLength() - 1;

   if ( last_index == -1 )
      return first_velocity;
   else
      return G4ThreeVector(   buffer[last_index][MOMENTUM_SLOT],
                              buffer[last_index][MOMENTUM_SLOT + 1],
                              buffer[last_index][MOMENTUM_SLOT + 2] ).mag();
}


void StepTracker::update_time_arclength( G4double time_to_add, G4double arclength_to_add) {

   // Right now assuming pure magnetic field, so velocity doesn't change.

   // If update_time_arclength() is called from context of position/momentum coordinates,
   // not position/velocity coordinates, it is the responsibility of the caller to
   // convert to position/velocity coordinates ( divide momentum coordinates
   // through by velocity).

   if ( ! within_AdvanceChordLimited ) {
      //cout << "Not within AdvancedChordLimited" << endl;
      return;
   }

   last_time_val_was_accepted = true;

   last_curve_length += arclength_to_add;
   last_time_length += time_to_add;


   G4int last_index = buffer.size() - 1;

   if (last_curve_length != buffer[last_index][ENDPOINT_BASE_INDEX + 1])
      cout << last_curve_length << ", " << buffer[last_index][ENDPOINT_BASE_INDEX + 1] << endl;
   if (last_time_length != buffer[last_index][ENDPOINT_BASE_INDEX + 0])
         cout << last_time_length << ", " << buffer[last_index][ENDPOINT_BASE_INDEX + 0] << endl;

   assert( last_curve_length == buffer[last_index][ENDPOINT_BASE_INDEX + 1] );
   assert( last_time_length == buffer[last_index][ENDPOINT_BASE_INDEX + 0] );



   /*
   G4int last_index = buffer.size() - 1;

   buffer[last_index][ENDPOINT_BASE_INDEX + 0] = last_time_length = buffer[last_index][0] + time_to_add;
   buffer[last_index][ENDPOINT_BASE_INDEX + 1] = last_curve_length = buffer[last_index][1] + arclength_to_add;
   */
}


