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

#ifndef MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYTIME_HH_
#define MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYTIME_HH_

#include "G4Mag_EqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ThreeVector.hh"
#include "G4CachedMagneticField.hh"
#include <vector>
#include "isTracking.hh"
#define NO_STATE_VARIABLES 12 // Used as upper bound for statically allocated array

using namespace std;

template <class BaseStepper>
class MagIntegratorStepper_byTime : public BaseStepper {
public:
   inline
   MagIntegratorStepper_byTime(G4Mag_EqRhs *EquationRhs,
                                 G4int numberOfVariables = 6,
                                 G4int numStateVariables = 12);
   virtual ~MagIntegratorStepper_byTime();

   inline void ComputeRightHandSide(const G4double yInput[], G4double dydx[]);



   inline void Stepper(const G4double yInput[],
                       const G4double dydx[],
                             G4double hstep,
                             G4double yOutput[],
                             G4double yError [] );

   inline G4double DistChord() const;

private:
   // Needed because Stepper() takes yInput[] and dydx[] as const input
   // but we have to do some scaling modifications in Stepper():
    // We want these to be statically allocated for faster access.
    // NO_STATE_VARIABLES is a functioning upper bound to num_variables.
   G4double yIn[NO_STATE_VARIABLES],
            dydx_copy[NO_STATE_VARIABLES];

   // Possible future use:
   //G4double cached_dydx[NO_STATE_VARIABLES],
   //         last_function_evaluation[NO_STATE_VARIABLES];

   G4double nextFunctionEvaluation[NO_STATE_VARIABLES];

   G4Mag_EqRhs  *m_fEq;
    int num_variables;
    int momentum_variables_index_offset;
};


template <class BaseStepper>
inline
void MagIntegratorStepper_byTime<BaseStepper>::ComputeRightHandSide(
                                             const G4double yInput[], G4double dydx[] ) {

   for (int i = 0; i < num_variables; i ++)
      yIn[i] = yInput[i];

   G4double current_relativistic_mass = m_fEq -> FMass();
   G4double rel_mass_inverse = 1. /current_relativistic_mass;
   for (int i = momentum_variables_index_offset; i < num_variables; i ++)
      yIn[i] *= rel_mass_inverse;

   BaseStepper::ComputeRightHandSide(yIn, dydx);

   for (int i = momentum_variables_index_offset; i < num_variables; i ++) {
                                       // may would usually want to make this from i = 0,
                                       // (but it really doesn't matter since we're only
                                       // using this for Nystrom Steppers.
      dydx[i] *= current_relativistic_mass;
   }
   // We always feed this template class a Mag_UsualEqRhs_IntegrateByTime as EquationRhs.
}

template <class BaseStepper>
inline
void MagIntegratorStepper_byTime<BaseStepper>::Stepper(const G4double yInput[],
      const G4double dydx[],
            G4double hstep,
            G4double yOutput[],
            G4double yError [] ) {

   // Have to copy because yInput and dydx are constant in the function signature...
   for (int i = 0; i < num_variables; i ++) {
      yIn[i] = yInput[i];
      dydx_copy[i] = dydx[i];
   }

   G4double current_relativistic_mass = m_fEq -> FMass();
   G4double rel_mass_inverse = 1. /current_relativistic_mass;
   // ...because now we have to convert to velocity coordinates:
   for (int i = momentum_variables_index_offset; i < num_variables; i ++) {
      yIn[i] *= rel_mass_inverse;
      dydx_copy[i] *= rel_mass_inverse;
   }

   /////// Getting velocity for purposes of rescaling the step length:
   // This is something which is specific to only one particle
   // so we don't use momentum_variables_index_offset or num_variables:
   G4double velocity = sqrt( yIn[3]*yIn[3] + yIn[4]*yIn[4] + yIn[5]*yIn[5] );

   // Within stepper call, convert hstep (which is in units of arclength)
   // to hstep / velocity (which is units of time).
   BaseStepper::Stepper( yIn, dydx_copy, hstep / velocity, yOutput, yError );

#ifdef TRACKING
   if ( BaseStepper::mTracker -> get_within_AdvanceChordLimited() ) {
      // Within AdvancedChordLimited, so we want to record.

      if ( BaseStepper::mTracker -> isArmed() ) { // If armed then store step result.

         // We need the next RHS function evaluation (for the right endpoint of
         // step interval). We need this because it currently is not stored by
         // the stepper (FSAL?)
         //
         // Since yOutput is currently in velocity coordinates, we just call
         // the BaseStepper::ComputeRightHandSide() method.
         BaseStepper::ComputeRightHandSide(yOutput, nextFunctionEvaluation);

         // previous function call should not be counted as part of total function calls:
         BaseStepper::mTracker -> add_to_num_other_function_calls_not_to_count(1);

         // Getting number of function calls used so far:
         const G4CachedMagneticField *myField = (G4CachedMagneticField*)
                           ( BaseStepper::GetEquationOfMotion() -> GetFieldObj() );
         G4int no_function_calls = myField -> GetCountCalls();

         // Storing all this with StepTracker:
         BaseStepper::mTracker -> RecordResultOfStepper(yIn, dydx_copy,
                                                        yOutput, nextFunctionEvaluation,
                                                        hstep, // Supposed to be in arclength units.
                                                        no_function_calls);

         BaseStepper::mTracker -> UnArmTracker();
      }
   }
#endif

   // Have to convert back to momentum coordinates.
   for (int i = momentum_variables_index_offset; i < num_variables; i ++) {
      yOutput[i] *= current_relativistic_mass;
      yError[i] *= current_relativistic_mass;      // Temp check to see if we have to scale last 3 coordinates of the error also.
   }
}

template <class BaseStepper>
inline
G4double MagIntegratorStepper_byTime<BaseStepper>::DistChord() const{

#ifdef NO_COUNT_FUNCTION_CALLS_FROM_DISTCHORD
   G4int no_function_calls_before_aux_stepper =
           (( G4CachedMagneticField* )( m_fEq -> GetFieldObj() ))
                                                  -> GetCountCalls();

#endif

   G4double dist_chord = BaseStepper::DistChord();

#ifdef NO_COUNT_FUNCTION_CALLS_FROM_DISTCHORD

   BaseStepper::mTracker -> add_to_num_other_function_calls_not_to_count(
          (( G4CachedMagneticField* )( m_fEq -> GetFieldObj() ))
             -> GetCountCalls() - no_function_calls_before_aux_stepper );
#endif

   return dist_chord;
}

template <class BaseStepper>
inline MagIntegratorStepper_byTime<BaseStepper>::MagIntegratorStepper_byTime(
                                                         G4Mag_EqRhs *EquationRhs,
                                                         G4int numberOfVariables,
                                                         G4int numStateVariables)
: BaseStepper( EquationRhs, numberOfVariables )
{
   num_variables = numberOfVariables;
   momentum_variables_index_offset = num_variables / 2;  // This assumes we have half position
                                                         // and half momentum variables
   m_fEq = EquationRhs;

   for (int i = 0; i < num_variables; i ++)
      yIn[i] = 0.;
}
template <class BaseStepper>
inline MagIntegratorStepper_byTime<BaseStepper>::~MagIntegratorStepper_byTime() {
}

#endif /* MAGNETICFIELD_INCLUDE_MAGINTEGRATORSTEPPERBYTIME_HH_ */

