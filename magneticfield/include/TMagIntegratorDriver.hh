#ifndef TMagInt_Driver_Def
#define TMagInt_Driver_Def

#include "G4Types.hh"
#include "G4FieldTrack.hh"
#include <iomanip>
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeometryTolerance.hh"

#ifndef G4NO_FIELD_STATISTICS
#define G4FLD_STATS  1
#endif

template
<class T_Stepper>
class TMagInt_Driver : public G4MagInt_Driver
{


    //  Stepsize can increase by no more than 5.0
    //           and decrease by no more than 1/10. = 0.1
    //
    static const G4double max_stepping_increase = 5.0;
    static const G4double max_stepping_decrease = 0.1;

    //  The (default) maximum number of steps is Base
    //  divided by the order of Stepper
    //
    static const G4int  fMaxStepBase = 250;  // Was 5000

    public:  // with description


    G4bool
        AccurateAdvance(G4FieldTrack& y_current,
                G4double     hstep,
                G4double     eps,
                G4double hinitial=0.0 )
        {
            // Runge-Kutta driver with adaptive stepsize control. Integrate starting
            // values at y_current over hstep x2 with accuracy eps. 
            // On output ystart is replaced by values at the end of the integration 
            // interval. RightHandSide is the right-hand side of ODE system. 
            // The source is similar to odeint routine from NRC p.721-722 .

            G4int nstp, i, no_warnings=0;
            G4double x, hnext, hdid, h;

#ifdef G4DEBUG_FIELD
            static G4int dbg=1;
            static G4int nStpPr=50;   // For debug printing of long integrations
            G4double ySubStepStart[G4FieldTrack::ncompSVEC];
            G4FieldTrack  yFldTrkStart(y_current);
#endif

            G4double y[G4FieldTrack::ncompSVEC], dydx[G4FieldTrack::ncompSVEC];
            G4double ystart[G4FieldTrack::ncompSVEC], yEnd[G4FieldTrack::ncompSVEC]; 
            G4double  x1, x2;
            G4bool succeeded = true, lastStepSucceeded;

            G4double startCurveLength;

            G4int  noFullIntegr=0, noSmallIntegr = 0 ;
            static G4ThreadLocal G4int  noGoodSteps =0 ;  // Bad = chord > curve-len 
            const  G4int  nvar= fNoVars;

            G4FieldTrack yStartFT(y_current);

            //  Ensure that hstep > 0
            //
            if( hstep <= 0.0 )
            { 
                if(hstep==0.0)
                {
                    std::ostringstream message;
                    message << "Proposed step is zero; hstep = " << hstep << " !";
                    G4Exception("AccurateAdvance()", 
                            "GeomField1001", JustWarning, message);
                    return succeeded; 
                }
                else
                { 
                    std::ostringstream message;
                    message << "Invalid run condition." << G4endl
                        << "Proposed step is negative; hstep = " << hstep << "." << G4endl
                        << "Requested step cannot be negative! Aborting event.";
                    G4Exception("AccurateAdvance()", 
                            "GeomField0003", EventMustBeAborted, message);
                    return false;
                }
            }

            y_current.DumpToArray( ystart );

            startCurveLength= y_current.GetCurveLength();
            x1= startCurveLength; 
            x2= x1 + hstep;

            if ( (hinitial > 0.0) && (hinitial < hstep)
                    && (hinitial > perMillion * hstep) )
            {
                h = hinitial;
            }
            else  //  Initial Step size "h" defaults to the full interval
            {
                h = hstep;
            }

            x = x1;

            for (i=0;i<nvar;i++)  { y[i] = ystart[i]; }

            G4bool lastStep= false;
            nstp=1;

            do
            {
                G4ThreeVector StartPos( y[0], y[1], y[2] );

#ifdef G4DEBUG_FIELD
                G4double xSubStepStart= x; 
                for (i=0;i<nvar;i++)  { ySubStepStart[i] = y[i]; }
                yFldTrkStart.LoadFromArray(y, fNoIntegrationVariables);
                yFldTrkStart.SetCurveLength(x);
#endif

                // Old method - inline call to Equation of Motion
                //   pIntStepper->RightHandSide( y, dydx );
                // New method allows to cache field, or state (eg momentum magnitude)
                pIntStepper->T_Stepper::RightHandSide( y, dydx );
                fNoTotalSteps++;

                // Perform the Integration
                //      
                if( h > fMinimumStep )
                { 
                    OneGoodStep(y,dydx,x,h,eps,hdid,hnext) ;
                    //--------------------------------------
                    lastStepSucceeded= (hdid == h);   
#ifdef G4DEBUG_FIELD
                    if (dbg>2) {
                        PrintStatus( ySubStepStart, xSubStepStart, y, x, h,  nstp); // Only
                    }
#endif
                }
                else
                {
                    G4FieldTrack yFldTrk( G4ThreeVector(0,0,0), 
                            G4ThreeVector(0,0,0), 0., 0., 0., 0. );
                    G4double dchord_step, dyerr, dyerr_len;   // What to do with these ?
                    yFldTrk.LoadFromArray(y, fNoIntegrationVariables); 
                    yFldTrk.SetCurveLength( x );

                    QuickAdvance( yFldTrk, dydx, h, dchord_step, dyerr_len ); 
                    //-----------------------------------------------------

                    yFldTrk.DumpToArray(y);    

#ifdef G4FLD_STATS
                    fNoSmallSteps++; 
                    if ( dyerr_len > fDyerr_max)  { fDyerr_max= dyerr_len; }
                    fDyerrPos_smTot += dyerr_len;
                    fSumH_sm += h;  // Length total for 'small' steps
                    if (nstp<=1)  { fNoInitialSmallSteps++; }
#endif
#ifdef G4DEBUG_FIELD
                    if (dbg>1)
                    {
                        if(fNoSmallSteps<2) { PrintStatus(ySubStepStart, x1, y, x, h, -nstp); }
                        G4cout << "Another sub-min step, no " << fNoSmallSteps 
                            << " of " << fNoTotalSteps << " this time " << nstp << G4endl; 
                        PrintStatus( ySubStepStart, x1, y, x, h,  nstp);   // Only this
                        G4cout << " dyerr= " << dyerr_len << " relative = " << dyerr_len / h 
                            << " epsilon= " << eps << " hstep= " << hstep 
                            << " h= " << h << " hmin= " << fMinimumStep << G4endl;
                    }
#endif        
                    if( h == 0.0 )
                    { 
                        G4Exception("AccurateAdvance()",
                                "GeomField0003", FatalException,
                                "Integration Step became Zero!"); 
                    }
                    dyerr = dyerr_len / h;
                    hdid= h;
                    x += hdid;

                    // Compute suggested new step
                    hnext= ComputeNewStepSize( dyerr/eps, h);

                    // .. hnext= ComputeNewStepSize_WithinLimits( dyerr/eps, h);
                    lastStepSucceeded= (dyerr<= eps);
                }

                if (lastStepSucceeded)  { noFullIntegr++; }
                else                    { noSmallIntegr++; }

                G4ThreeVector EndPos( y[0], y[1], y[2] );

#ifdef  G4DEBUG_FIELD
                if( (dbg>0) && (dbg<=2) && (nstp>nStpPr))
                {
                    if( nstp==nStpPr )  { G4cout << "***** Many steps ****" << G4endl; }
                    G4cout << "MagIntDrv: " ; 
                    G4cout << "hdid="  << std::setw(12) << hdid  << " "
                        << "hnext=" << std::setw(12) << hnext << " " 
                        << "hstep=" << std::setw(12) << hstep << " (requested) " 
                        << G4endl;
                    PrintStatus( ystart, x1, y, x, h, (nstp==nStpPr) ? -nstp: nstp); 
                }
#endif

                // Check the endpoint
                G4double endPointDist= (EndPos-StartPos).mag(); 
                if ( endPointDist >= hdid*(1.+perMillion) )
                {
                    fNoBadSteps++;

                    // Issue a warning only for gross differences -
                    // we understand how small difference occur.
                    if ( endPointDist >= hdid*(1.+perThousand) )
                    { 
#ifdef G4DEBUG_FIELD
                        if (dbg)
                        {
                            WarnEndPointTooFar ( endPointDist, hdid, eps, dbg ); 
                            G4cerr << "  Total steps:  bad " << fNoBadSteps
                                << " good " << noGoodSteps << " current h= " << hdid
                                << G4endl;
                            PrintStatus( ystart, x1, y, x, hstep, no_warnings?nstp:-nstp);  
                        }
#endif
                        no_warnings++;
                    }
                }
                else
                {
                    noGoodSteps ++;
                } 
                // #endif

                //  Avoid numerous small last steps
                if( (h < eps * hstep) || (h < fSmallestFraction * startCurveLength) )
                {
                    // No more integration -- the next step will not happen
                    lastStep = true;  
                }
                else
                {
                    // Check the proposed next stepsize
                    if(std::fabs(hnext) <= Hmin())
                    {
#ifdef  G4DEBUG_FIELD
                        // If simply a very small interval is being integrated, do not warn
                        if( (x < x2 * (1-eps) ) &&        // The last step can be small: OK
                                (std::fabs(hstep) > Hmin()) ) // and if we are asked, it's OK
                        {
                            if(dbg>0)
                            {
                                WarnSmallStepSize( hnext, hstep, h, x-x1, nstp );  
                                PrintStatus( ystart, x1, y, x, hstep, no_warnings?nstp:-nstp);
                            }
                            no_warnings++;
                        }
#endif
                        // Make sure that the next step is at least Hmin.
                        h = Hmin();
                    }
                    else
                    {
                        h = hnext;
                    }

                    //  Ensure that the next step does not overshoot
                    if ( x+h > x2 )
                    {                // When stepsize overshoots, decrease it!
                        h = x2 - x ;   // Must cope with difficult rounding-error
                    }                // issues if hstep << x2

                    if ( h == 0.0 )
                    {
                        // Cannot progress - accept this as last step - by default
                        lastStep = true;
#ifdef G4DEBUG_FIELD
                        if (dbg>2)
                        {
                            int prec= G4cout.precision(12); 
                            G4cout << "Warning: TMagIntegratorDriver::AccurateAdvance"
                                << G4endl
                                << "  Integration step 'h' became "
                                << h << " due to roundoff. " << G4endl
                                << " Calculated as difference of x2= "<< x2 << " and x=" << x
                                << "  Forcing termination of advance." << G4endl;
                            G4cout.precision(prec);
                        }          
#endif
                    }
                }
            } while ( ((nstp++)<=fMaxNoSteps) && (x < x2) && (!lastStep) );
            // Have we reached the end ?
            // --> a better test might be x-x2 > an_epsilon

            succeeded=  (x>=x2);  // If it was a "forced" last step

            for (i=0;i<nvar;i++)  { yEnd[i] = y[i]; }

            // Put back the values.
            y_current.LoadFromArray( yEnd, fNoIntegrationVariables );
            y_current.SetCurveLength( x );

            if(nstp > fMaxNoSteps)
            {
                no_warnings++;
                succeeded = false;
#ifdef G4DEBUG_FIELD
                if (dbg)
                {
                    WarnTooManyStep( x1, x2, x );  //  Issue WARNING
                    PrintStatus( yEnd, x1, y, x, hstep, -nstp);
                }
#endif
            }

#ifdef G4DEBUG_FIELD
            if( dbg && no_warnings )
            {
                G4cerr << "TMagIntegratorDriver exit status: no-steps " << nstp <<G4endl;
                PrintStatus( yEnd, x1, y, x, hstep, nstp);
            }
#endif

            return succeeded;
        }  // end of AccurateAdvance ...........................

    // ---------------------------------------------------------

    // QuickAdvance just tries one Step - it does not ensure accuracy
    //
    G4bool  QuickAdvance(       
            G4FieldTrack& y_posvel,         // INOUT
            const G4double     dydx[],  
            G4double     hstep,       // In
            G4double&    dchord_step,
            G4double&    dyerr_pos_sq,
            G4double&    dyerr_mom_rel_sq )  
    {
        G4Exception("QuickAdvance()", "GeomField0001",
                FatalException, "Not yet implemented."); 

        // Use the parameters of this method, to please compiler
        dchord_step = dyerr_pos_sq = hstep * hstep * dydx[0]; 
        dyerr_mom_rel_sq = y_posvel.GetPosition().mag2();
        return true;
    }

    //----------------------------------------------------------------------

    G4bool  QuickAdvance(       
            G4FieldTrack& y_posvel,         // INOUT
            const G4double     dydx[],  
            G4double     hstep,       // In
            G4double&    dchord_step,
            G4double&    dyerr )
    {
        G4double dyerr_pos_sq, dyerr_mom_rel_sq;  
        G4double yerr_vec[G4FieldTrack::ncompSVEC],
                 yarrin[G4FieldTrack::ncompSVEC], yarrout[G4FieldTrack::ncompSVEC]; 
        G4double s_start;
        G4double dyerr_mom_sq, vel_mag_sq, inv_vel_mag_sq;

        static G4ThreadLocal G4int no_call=0; 
        no_call ++; 

        // Move data into array
        y_posvel.DumpToArray( yarrin );      //  yarrin  <== y_posvel 
        s_start = y_posvel.GetCurveLength();

        // Do an Integration Step
        pIntStepper-> T_Stepper::Stepper(yarrin, dydx, hstep, yarrout, yerr_vec) ; 
        //            *******

        // Estimate curve-chord distance
        dchord_step= pIntStepper-> T_Stepper::DistChord();
        //                         *********

        // Put back the values.  yarrout ==> y_posvel
        y_posvel.LoadFromArray( yarrout, fNoIntegrationVariables );
        y_posvel.SetCurveLength( s_start + hstep );

#ifdef  G4DEBUG_FIELD
        if(fVerboseLevel>2)
        {
            G4cout << "G4MagIntDrv: Quick Advance" << G4endl;
            PrintStatus( yarrin, s_start, yarrout, s_start+hstep, hstep,  1); 
        }
#endif

        // A single measure of the error   
        //      TO-DO :  account for  energy,  spin, ... ? 
        vel_mag_sq   = ( sqr(yarrout[3])+sqr(yarrout[4])+sqr(yarrout[5]) );
        inv_vel_mag_sq = 1.0 / vel_mag_sq; 
        dyerr_pos_sq = ( sqr(yerr_vec[0])+sqr(yerr_vec[1])+sqr(yerr_vec[2]));
        dyerr_mom_sq = ( sqr(yerr_vec[3])+sqr(yerr_vec[4])+sqr(yerr_vec[5]));
        dyerr_mom_rel_sq = dyerr_mom_sq * inv_vel_mag_sq;

        // Calculate also the change in the momentum squared also ???
        // G4double veloc_square = y_posvel.GetVelocity().mag2();
        // ...

#ifdef RETURN_A_NEW_STEP_LENGTH
        // The following step cannot be done here because "eps" is not known.
        dyerr_len = std::sqrt( dyerr_len_sq ); 
        dyerr_len_sq /= eps ;

        // Look at the velocity deviation ?
        //  sqr(yerr_vec[3])+sqr(yerr_vec[4])+sqr(yerr_vec[5]));

        // Set suggested new step
        hstep= ComputeNewStepSize( dyerr_len, hstep);
#endif

        if( dyerr_pos_sq > ( dyerr_mom_rel_sq * sqr(hstep) ) )
        {
            dyerr = std::sqrt(dyerr_pos_sq);
        }
        else
        {
            // Scale it to the current step size - for now
            dyerr = std::sqrt(dyerr_mom_rel_sq) * hstep;
        }

        return true;
    }

    // --------------------------------------------------------------------------

    //  Constructor
    //

    TMagInt_Driver( G4double                hminimum, 
            T_Stepper *pStepper,
            G4int                   numComponents=6,
            G4int                   statisticsVerbose=1)
        : fSmallestFraction( 1.0e-12 ), 
        fNoIntegrationVariables(numComponents), 
        fMinNoVars(12), 
        fNoVars( std::max( fNoIntegrationVariables, fMinNoVars )),
        fStatisticsVerboseLevel(statisticsVerbose),
        fNoTotalSteps(0),  fNoBadSteps(0), fNoSmallSteps(0),
        fNoInitialSmallSteps(0), 
        fDyerr_max(0.0), fDyerr_mx2(0.0), 
        fDyerrPos_smTot(0.0), fDyerrPos_lgTot(0.0), fDyerrVel_lgTot(0.0), 
        fSumH_sm(0.0), fSumH_lg(0.0),
        fVerboseLevel(0), 
        G4MagInt_Driver(hminimum, pStepper, numComponents, statisticsVerbose)
    {  
        // In order to accomodate "Laboratory Time", which is [7], fMinNoVars=8
        // is required. For proper time of flight and spin,  fMinNoVars must be 12

        RenewStepperAndAdjust( pStepper );
        fMinimumStep= hminimum;
        fMaxNoSteps = fMaxStepBase / pIntStepper->T_Stepper::IntegratorOrder();
#ifdef G4DEBUG_FIELD
        fVerboseLevel=2;
#endif

        if( (fVerboseLevel > 0) || (fStatisticsVerboseLevel > 1) )
        {
            G4cout << "MagIntDriver version: Accur-Adv: "
                << "invE_nS, QuickAdv-2sqrt with Statistics "
#ifdef G4FLD_STATS
                << " enabled "
#else
                << " disabled "
#endif
                << G4endl;
        }
    }

    // ---------------------------------------------------------
    //  Destructor
    //
    ~TMagInt_Driver()
    { 
        if( fStatisticsVerboseLevel > 1 )
        {
            PrintStatisticsReport();
        }
    }

    // To add much printing for debugging purposes, uncomment the following
    // and set verbose level to 1 or higher value !
    // #define  G4DEBUG_FIELD 1    

    // ---------------------------------------------------------



    inline
        G4double GetHmin() const
        {
            return fMinimumStep;
        } 

    inline
        G4double Hmin() const
        {
            return fMinimumStep;
        }

    inline
        G4double GetSafety() const
        {
            return safety;
        }

    inline
        G4double GetPshrnk() const
        {
            return pshrnk;
        } 

    inline
        G4double GetPgrow() const
        {
            return pgrow;
        }

    inline
        G4double GetErrcon() const
        {
            return errcon;
        }

    inline
        void SetHmin(G4double newval)
        {
            fMinimumStep = newval;
        } 

    inline
        G4double ComputeAndSetErrcon()
        {
            errcon = std::pow(max_stepping_increase/GetSafety(),1.0/GetPgrow());
            return errcon;
        } 

    inline
        void ReSetParameters(G4double new_safety=0.9)
        {
            safety = new_safety;
            pshrnk = -1.0 / pIntStepper->T_Stepper::IntegratorOrder();
            pgrow  = -1.0 / (1.0 + pIntStepper->T_Stepper::IntegratorOrder());
            ComputeAndSetErrcon();
        }

    inline
        void SetSafety(G4double val)
        { 
            safety=val;
            ComputeAndSetErrcon();
        }

    inline
        void SetPgrow(G4double  val)
        { 
            pgrow=val;
            ComputeAndSetErrcon(); 
        }

    inline
        void SetErrcon(G4double val)
        { 
            errcon=val;
        }

    inline
        void RenewStepperAndAdjust(T_Stepper *pItsStepper)
        {  
            pIntStepper = pItsStepper; 
            ReSetParameters();
        }

    inline
        const T_Stepper* GetStepper() const
        {
            return pIntStepper;
        }

    inline
        T_Stepper* GetStepper() 
        {
            return pIntStepper;
        }

    inline
        G4int GetMaxNoSteps() const
        {
            return fMaxNoSteps;
        }

    inline
        void SetMaxNoSteps(G4int val)
        {
            fMaxNoSteps= val;
        }

    inline
        void GetDerivatives(const G4FieldTrack &y_curr, // const, INput
                G4double     dydx[])  // OUTput
        { 
            G4double  tmpValArr[G4FieldTrack::ncompSVEC];
            y_curr.DumpToArray( tmpValArr  );
            pIntStepper -> T_Stepper::RightHandSide( tmpValArr , dydx );
        }

    inline
        G4double GetVerboseLevel() const
        {
            return fVerboseLevel;
        } 

    inline 
        void SetVerboseLevel(G4int newLevel)
        {
            fVerboseLevel= newLevel;
        }

    inline
        G4double GetSmallestFraction() const
        {
            return fSmallestFraction; 
        } 

    void
        OneGoodStep(      G4double y[],        // InOut
                const G4double dydx[],
                G4double& x,         // InOut
                G4double htry,
                G4double eps_rel_max,
                G4double& hdid,      // Out
                G4double& hnext )    // Out

        // Driver for one Runge-Kutta Step with monitoring of local truncation error
        // to ensure accuracy and adjust stepsize. Input are dependent variable
        // array y[0,...,5] and its derivative dydx[0,...,5] at the
        // starting value of the independent variable x . Also input are stepsize
        // to be attempted htry, and the required accuracy eps. On output y and x
        // are replaced by their new values, hdid is the stepsize that was actually
        // accomplished, and hnext is the estimated next stepsize. 
        // This is similar to the function rkqs from the book:
        // Numerical Recipes in C: The Art of Scientific Computing (NRC), Second
        // Edition, by William H. Press, Saul A. Teukolsky, William T.
        // Vetterling, and Brian P. Flannery (Cambridge University Press 1992),
        // 16.2 Adaptive StepSize Control for Runge-Kutta, p. 719

        {
            G4double errmax_sq;
            G4double h, htemp, xnew ;

            G4double yerr[G4FieldTrack::ncompSVEC], ytemp[G4FieldTrack::ncompSVEC];

            h = htry ; // Set stepsize to the initial trial value

            G4double inv_eps_vel_sq = 1.0 / (eps_rel_max*eps_rel_max);

            G4double errpos_sq=0.0;    // square of displacement error
            G4double errvel_sq=0.0;    // square of momentum vector difference
            G4double errspin_sq=0.0;   // square of spin vector difference

            G4int iter;

            static G4ThreadLocal G4int tot_no_trials=0; 
            const G4int max_trials=100; 

            G4ThreeVector Spin(y[9],y[10],y[11]);
            G4bool     hasSpin= (Spin.mag2() > 0.0); 

            for (iter=0; iter<max_trials ;iter++)
            {
                tot_no_trials++;
                pIntStepper-> T_Stepper::Stepper(y,dydx,h,ytemp,yerr); 
                //            *******
                G4double eps_pos = eps_rel_max * std::max(h, fMinimumStep); 
                G4double inv_eps_pos_sq = 1.0 / (eps_pos*eps_pos); 

                // Evaluate accuracy
                //
                errpos_sq =  sqr(yerr[0]) + sqr(yerr[1]) + sqr(yerr[2]) ;
                errpos_sq *= inv_eps_pos_sq; // Scale relative to required tolerance

                // Accuracy for momentum
                errvel_sq =  (sqr(yerr[3]) + sqr(yerr[4]) + sqr(yerr[5]) )
                    / (sqr(y[3]) + sqr(y[4]) + sqr(y[5]) );
                errvel_sq *= inv_eps_vel_sq;
                errmax_sq = std::max( errpos_sq, errvel_sq ); // Square of maximum error

                if( hasSpin )
                { 
                    // Accuracy for spin
                    errspin_sq =  ( sqr(yerr[9]) + sqr(yerr[10]) + sqr(yerr[11]) )
                        /  ( sqr(y[9]) + sqr(y[10]) + sqr(y[11]) );
                    errspin_sq *= inv_eps_vel_sq;
                    errmax_sq = std::max( errmax_sq, errspin_sq ); 
                }

                if ( errmax_sq <= 1.0 )  { break; } // Step succeeded. 

                // Step failed; compute the size of retrial Step.
                htemp = GetSafety()*h* std::pow( errmax_sq, 0.5*GetPshrnk() );

                if (htemp >= 0.1*h)  { h = htemp; }  // Truncation error too large,
                else  { h = 0.1*h; }                 // reduce stepsize, but no more
                // than a factor of 10
                xnew = x + h;
                if(xnew == x)
                {
                    G4cerr << "G4MagIntegratorDriver::OneGoodStep:" << G4endl
                        << "  Stepsize underflow in Stepper " << G4endl ;
                    G4cerr << "  Step's start x=" << x << " and end x= " << xnew 
                        << " are equal !! " << G4endl
                        <<"  Due to step-size= " << h 
                        << " . Note that input step was " << htry << G4endl;
                    break;
                }
            }

#ifdef G4FLD_STATS
            // Sum of squares of position error // and momentum dir (underestimated)
            fSumH_lg += h; 
            fDyerrPos_lgTot += errpos_sq;
            fDyerrVel_lgTot += errvel_sq * h * h; 
#endif

            // Compute size of next Step
            if (errmax_sq > errcon*errcon)
            { 
                hnext = GetSafety()*h*std::pow(errmax_sq, 0.5*GetPgrow());
            }
            else
            {
                hnext = max_stepping_increase*h ; // No more than a factor of 5 increase
            }
            x += (hdid = h);

            for(G4int k=0;k<fNoIntegrationVariables;k++) { y[k] = ytemp[k]; }

            return;
        }   // end of  OneGoodStep .............................

    //----------------------------------------------------------------------

    //  This method computes new step sizes - but does not limit changes to
    //   within  certain factors
    // 
    G4double 
        ComputeNewStepSize( 
                G4double  errMaxNorm,    // max error  (normalised)
                G4double  hstepCurrent)  // current step size
        {
            G4double hnew;

            // Compute size of next Step for a failed step
            if(errMaxNorm > 1.0 )
            {
                // Step failed; compute the size of retrial Step.
                hnew = GetSafety()*hstepCurrent*std::pow(errMaxNorm,GetPshrnk()) ;
            } else if(errMaxNorm > 0.0 ) {
                // Compute size of next Step for a successful step
                hnew = GetSafety()*hstepCurrent*std::pow(errMaxNorm,GetPgrow()) ;
            } else {
                // if error estimate is zero (possible) or negative (dubious)
                hnew = max_stepping_increase * hstepCurrent; 
            }

            return hnew;
        }

    // ---------------------------------------------------------------------------

    // This method computes new step sizes limiting changes within certain factors
    // 
    // It shares its logic with AccurateAdvance.
    // They are kept separate currently for optimisation.
    //
    G4double 
        ComputeNewStepSize_WithinLimits( 
                G4double  errMaxNorm,    // max error  (normalised)
                G4double  hstepCurrent)  // current step size
        {
            G4double hnew;

            // Compute size of next Step for a failed step
            if (errMaxNorm > 1.0 )
            {
                // Step failed; compute the size of retrial Step.
                hnew = GetSafety()*hstepCurrent*std::pow(errMaxNorm,GetPshrnk()) ;

                if (hnew < max_stepping_decrease*hstepCurrent)
                {
                    hnew = max_stepping_decrease*hstepCurrent ;
                    // reduce stepsize, but no more
                    // than this factor (value= 1/10)
                }
            }
            else
            {
                // Compute size of next Step for a successful step
                if (errMaxNorm > errcon)
                { hnew = GetSafety()*hstepCurrent*std::pow(errMaxNorm,GetPgrow()); }
                else  // No more than a factor of 5 increase
                { hnew = max_stepping_increase * hstepCurrent; }
            }
            return hnew;
        }

    // ---------------------------------------------------------------------------



    public:  // without description


    void SetSmallestFraction(G4double newFraction)
    {
        if( (newFraction > 1.e-16) && (newFraction < 1e-8) )
        {
            fSmallestFraction= newFraction;
        }
        else
        { 
            G4cerr << "Warning: SmallestFraction not changed. " << G4endl
                << "  Proposed value was " << newFraction << G4endl
                << "  Value must be between 1.e-8 and 1.e-16" << G4endl;
        }
    }

    protected:  // without description
    void
        WarnSmallStepSize( G4double hnext, G4double hstep, 
                G4double h, G4double xDone,
                G4int nstp)
        {
            static G4ThreadLocal G4int noWarningsIssued =0;
            const  G4int maxNoWarnings =  10;   // Number of verbose warnings
            std::ostringstream message;
            if( (noWarningsIssued < maxNoWarnings) || fVerboseLevel > 10 )
            {
                message << "The stepsize for the next iteration, " << hnext
                    << ", is too small - in Step number " << nstp << "." << G4endl
                    << "The minimum for the driver is " << Hmin()  << G4endl
                    << "Requested integr. length was " << hstep << " ." << G4endl
                    << "The size of this sub-step was " << h     << " ." << G4endl
                    << "The integrations has already gone " << xDone;
            }
            else
            {
                message << "Too small 'next' step " << hnext
                    << ", step-no: " << nstp << G4endl
                    << ", this sub-step: " << h     
                    << ",  req_tot_len: " << hstep 
                    << ", done: " << xDone << ", min: " << Hmin();
            }
            G4Exception("WarnSmallStepSize()", "GeomField1001",
                    JustWarning, message);
            noWarningsIssued++;
        }

    // ---------------------------------------------------------
    void
        WarnTooManyStep( G4double x1start, 
                G4double x2end, 
                G4double xCurrent)
        {
            std::ostringstream message;
            message << "The number of steps used in the Integration driver"
                << " (Runge-Kutta) is too many." << G4endl
                << "Integration of the interval was not completed !" << G4endl
                << "Only a " << (xCurrent-x1start)*100/(x2end-x1start)
                << " % fraction of it was done.";
            G4Exception("WarnTooManyStep()", "GeomField1001",
                    JustWarning, message);
        }

    // ---------------------------------------------------------
    void
        WarnEndPointTooFar (G4double endPointDist, 
                G4double   h , 
                G4double  eps,
                G4int     dbg)
        {
            static G4ThreadLocal G4double maxRelError=0.0;
            G4bool isNewMax, prNewMax;

            isNewMax = endPointDist > (1.0 + maxRelError) * h;
            prNewMax = endPointDist > (1.0 + 1.05 * maxRelError) * h;
            if( isNewMax ) { maxRelError= endPointDist / h - 1.0; }

            if( dbg && (h > G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()) 
                    && ( (dbg>1) || prNewMax || (endPointDist >= h*(1.+eps) ) ) )
            { 
                static G4ThreadLocal G4int noWarnings = 0;
                std::ostringstream message;
                if( (noWarnings ++ < 10) || (dbg>2) )
                {
                    message << "The integration produced an end-point which " << G4endl
                        << "is further from the start-point than the curve length."
                        << G4endl;
                }
                message << "  Distance of endpoints = " << endPointDist
                    << ", curve length = " << h << G4endl
                    << "  Difference (curveLen-endpDist)= " << (h - endPointDist)
                    << ", relative = " << (h-endPointDist) / h 
                    << ", epsilon =  " << eps;
                G4Exception("WarnEndPointTooFar()", "GeomField1001",
                        JustWarning, message);
            }
        }

    // ---------------------------------------------------------
    // ---------------------------------------------------------------------------

    void PrintStatus( const G4double*   StartArr,  
            G4double          xstart,
            const G4double*   CurrentArr, 
            G4double          xcurrent,
            G4double          requestStep, 
            G4int             subStepNo)
        // Potentially add as arguments:  
        //                                 <dydx>           - as Initial Force
        //                                 stepTaken(hdid)  - last step taken
        //                                 nextStep (hnext) - proposal for size
    {
        G4FieldTrack  StartFT(G4ThreeVector(0,0,0),
                G4ThreeVector(0,0,0), 0., 0., 0., 0. );
        G4FieldTrack  CurrentFT (StartFT);

        StartFT.LoadFromArray( StartArr, fNoIntegrationVariables); 
        StartFT.SetCurveLength( xstart);
        CurrentFT.LoadFromArray( CurrentArr, fNoIntegrationVariables); 
        CurrentFT.SetCurveLength( xcurrent );

        PrintStatus(StartFT, CurrentFT, requestStep, subStepNo ); 
    }

    // ---------------------------------------------------------------------------

    void PrintStatus(
            const G4FieldTrack&  StartFT,
            const G4FieldTrack&  CurrentFT, 
            G4double             requestStep, 
            G4int                subStepNo)
    {
        G4int verboseLevel= fVerboseLevel;
        static G4ThreadLocal G4int noPrecision= 5;
        G4int oldPrec= G4cout.precision(noPrecision);
        // G4cout.setf(ios_base::fixed,ios_base::floatfield);

        const G4ThreeVector StartPosition=       StartFT.GetPosition();
        const G4ThreeVector StartUnitVelocity=   StartFT.GetMomentumDir();
        const G4ThreeVector CurrentPosition=     CurrentFT.GetPosition();
        const G4ThreeVector CurrentUnitVelocity= CurrentFT.GetMomentumDir();

        G4double  DotStartCurrentVeloc= StartUnitVelocity.dot(CurrentUnitVelocity);

        G4double step_len= CurrentFT.GetCurveLength() - StartFT.GetCurveLength();
        G4double subStepSize = step_len;

        if( (subStepNo <= 1) || (verboseLevel > 3) )
        {
            subStepNo = - subStepNo;        // To allow printing banner

            G4cout << std::setw( 6)  << " " << std::setw( 25)
                << " TMagInt_Driver: Current Position  and  Direction" << " "
                << G4endl; 
            G4cout << std::setw( 5) << "Step#" << " "
                << std::setw( 7) << "s-curve" << " "
                << std::setw( 9) << "X(mm)" << " "
                << std::setw( 9) << "Y(mm)" << " "  
                << std::setw( 9) << "Z(mm)" << " "
                << std::setw( 8) << " N_x " << " "
                << std::setw( 8) << " N_y " << " "
                << std::setw( 8) << " N_z " << " "
                << std::setw( 8) << " N^2-1 " << " "
                << std::setw(10) << " N(0).N " << " "
                << std::setw( 7) << "KinEner " << " "
                << std::setw(12) << "Track-l" << " "   // Add the Sub-step ??
                << std::setw(12) << "Step-len" << " " 
                << std::setw(12) << "Step-len" << " " 
                << std::setw( 9) << "ReqStep" << " "  
                << G4endl;
        }

        if( (subStepNo <= 0) )
        {
            PrintStat_Aux( StartFT,  requestStep, 0., 
                    0,        0.0,         1.0);
            //*************
        }

        if( verboseLevel <= 3 )
        {
            G4cout.precision(noPrecision);
            PrintStat_Aux( CurrentFT, requestStep, step_len, 
                    subStepNo, subStepSize, DotStartCurrentVeloc );
            //*************
        }

        else // if( verboseLevel > 3 )
        {
            //  Multi-line output

            // G4cout << "Current  Position is " << CurrentPosition << G4endl 
            //    << " and UnitVelocity is " << CurrentUnitVelocity << G4endl;
            // G4cout << "Step taken was " << step_len  
            //    << " out of PhysicalStep= " <<  requestStep << G4endl;
            // G4cout << "Final safety is: " << safety << G4endl;
            // G4cout << "Chord length = " << (CurrentPosition-StartPosition).mag()
            //        << G4endl << G4endl; 
        }
        G4cout.precision(oldPrec);
    }

    // ---------------------------------------------------------------------------

    void PrintStat_Aux(
            const G4FieldTrack&  aFieldTrack,
            G4double             requestStep, 
            G4double             step_len,
            G4int                subStepNo,
            G4double             subStepSize,
            G4double             dotVeloc_StartCurr)
    {
        const G4ThreeVector Position=      aFieldTrack.GetPosition();
        const G4ThreeVector UnitVelocity=  aFieldTrack.GetMomentumDir();

        if( subStepNo >= 0)
        {
            G4cout << std::setw( 5) << subStepNo << " ";
        }
        else
        {
            G4cout << std::setw( 5) << "Start" << " ";
        }
        G4double curveLen= aFieldTrack.GetCurveLength();
        G4cout << std::setw( 7) << curveLen;
        G4cout << std::setw( 9) << Position.x() << " "
            << std::setw( 9) << Position.y() << " "
            << std::setw( 9) << Position.z() << " "
            << std::setw( 8) << UnitVelocity.x() << " "
            << std::setw( 8) << UnitVelocity.y() << " "
            << std::setw( 8) << UnitVelocity.z() << " ";
        G4int oldprec= G4cout.precision(3);
        G4cout << std::setw( 8) << UnitVelocity.mag2()-1.0 << " ";
        G4cout.precision(6);
        G4cout << std::setw(10) << dotVeloc_StartCurr << " ";
        G4cout.precision(oldprec);
        G4cout << std::setw( 7) << aFieldTrack.GetKineticEnergy();
        G4cout << std::setw(12) << step_len << " ";

        static G4ThreadLocal G4double oldCurveLength= 0.0;
        static G4ThreadLocal G4double oldSubStepLength= 0.0;
        static G4ThreadLocal G4int oldSubStepNo= -1;

        G4double subStep_len=0.0;
        if( curveLen > oldCurveLength )
        {
            subStep_len= curveLen - oldCurveLength;
        }
        else if (subStepNo == oldSubStepNo)
        {
            subStep_len= oldSubStepLength;
        }
        oldCurveLength= curveLen;
        oldSubStepLength= subStep_len;

        G4cout << std::setw(12) << subStep_len << " "; 
        G4cout << std::setw(12) << subStepSize << " "; 
        if( requestStep != -1.0 )
        {
            G4cout << std::setw( 9) << requestStep << " ";
        }
        else
        {
            G4cout << std::setw( 9) << " InitialStep " << " ";
        }
        G4cout << G4endl;
    }

    // ---------------------------------------------------------------------------

    void PrintStatisticsReport()
    {
        G4int noPrecBig= 6;
        G4int oldPrec= G4cout.precision(noPrecBig);

        G4cout << "TMagInt_Driver Statistics of steps undertaken. " << G4endl;
        G4cout << "TMagInt_Driver: Number of Steps: "
            << " Total= " <<  fNoTotalSteps
            << " Bad= "   <<  fNoBadSteps 
            << " Small= " <<  fNoSmallSteps 
            << " Non-initial small= " << (fNoSmallSteps-fNoInitialSmallSteps)
            << G4endl;

#ifdef G4FLD_STATS
        G4cout << "MID dyerr: " 
            << " maximum= " << fDyerr_max 
            << " Sum small= " << fDyerrPos_smTot 
            << " std::sqrt(Sum large^2): pos= " << std::sqrt(fDyerrPos_lgTot)
            << " vel= " << std::sqrt( fDyerrVel_lgTot )
            << " Total h-distance: small= " << fSumH_sm 
            << " large= " << fSumH_lg
            << G4endl;

#if 0
        G4int noPrecSmall=4; 
        // Single line precis of statistics ... optional
        G4cout.precision(noPrecSmall);
        G4cout << "MIDnums: " << fMinimumStep
            << "   " << fNoTotalSteps 
            << "  "  <<  fNoSmallSteps
            << "  "  << fNoSmallSteps-fNoInitialSmallSteps
            << "  "  << fNoBadSteps         
            << "   " << fDyerr_max
            << "   " << fDyerr_mx2 
            << "   " << fDyerrPos_smTot 
            << "   " << fSumH_sm
            << "   " << fDyerrPos_lgTot
            << "   " << fDyerrVel_lgTot
            << "   " << fSumH_lg
            << G4endl;
#endif 
#endif 

        G4cout.precision(oldPrec);
    }

    // ---------------------------------------------------------------------------


#ifdef QUICK_ADV_TWO
    G4bool  QuickAdvance(       
            G4double     yarrin[],    // In
            const G4double     dydx[],  
            G4double     hstep,       // In
            G4double     yarrout[],
            G4double&    dchord_step,
            G4double&    dyerr )      // In length
    {
        G4Exception("QuickAdvance()", "GeomField0001",
                FatalException, "Not yet implemented.");
        dyerr = dchord_step = hstep * yarrin[0] * dydx[0];
        yarrout[0]= yarrin[0];
    }
#endif

    private:

    TMagInt_Driver(const TMagInt_Driver&);
    TMagInt_Driver& operator=(const TMagInt_Driver&);
    // Private copy constructor and assignment operator.

    private:

    // ---------------------------------------------------------------
    //  INVARIANTS 

    G4double  fMinimumStep;
    // Minimum Step allowed in a Step (in absolute units)
    G4double  fSmallestFraction;      //   Expected range 1e-12 to 5e-15;  
    // Smallest fraction of (existing) curve length - in relative units
    //  below this fraction the current step will be the last 

    const G4int  fNoIntegrationVariables;  // Number of Variables in integration
    const G4int  fMinNoVars;               // Minimum number for FieldTrack
    const G4int  fNoVars;                  // Full number of variable

    G4int   fMaxNoSteps;

    G4double safety;
    G4double pshrnk;   //  exponent for shrinking
    G4double pgrow;    //  exponent for growth
    G4double errcon;
    // Parameters used to grow and shrink trial stepsize.

    G4int    fStatisticsVerboseLevel;

    // ---------------------------------------------------------------
    // DEPENDENT Objects
    T_Stepper *pIntStepper;

    // ---------------------------------------------------------------
    //  STATE

    G4int  fNoTotalSteps, fNoBadSteps, fNoSmallSteps, fNoInitialSmallSteps; 
    G4double fDyerr_max, fDyerr_mx2;
    G4double fDyerrPos_smTot, fDyerrPos_lgTot, fDyerrVel_lgTot; 
    G4double fSumH_sm, fSumH_lg; 
    // Step Statistics 

    G4int  fVerboseLevel;   // Verbosity level for printing (debug, ..)
    // Could be varied during tracking - to help identify issues

};

#endif /* TMagInt_Driver_Def */
