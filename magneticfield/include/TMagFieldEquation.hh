

#include "G4ChargeState.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "sqrt.h"

#include <blaze/math/StaticVector.h>
#include <blaze/math/DenseSubvector.h>


using namespace blaze; 

typedef StaticVector<G4double, 3UL> Blaze3DVec;
typedef StaticVector<G4double, 8UL> Blaze8DVec;

template 
<class Field_t>
class TMagFieldEquation : public G4Mag_UsualEqRhs
{
    public:

        typedef Field_t T_Field;

        TMagFieldEquation(T_Field* f)
            : G4Mag_UsualEqRhs(f)
        {
            itsField = f;
        }

        virtual ~TMagFieldEquation(){;}

        __attribute__((always_inline)) 
            void GetFieldValue(const G4double Point[4],
                    G4double Field[]) const
            {
                itsField->T_Field::GetFieldValue(Point, Field);
            }

        __attribute__((always_inline)) 
            void RightHandSide(const G4double y[], G4double dydx[] )
            //	const
            {
                G4double Field[G4maximum_number_of_field_components]; 
                G4double  PositionAndTime[4];
                PositionAndTime[0] = y[0];
                PositionAndTime[1] = y[1];
                PositionAndTime[2] = y[2];
                PositionAndTime[3] = y[7];   
                GetFieldValue(PositionAndTime, Field) ;

                Blaze3DVec yv(y[3],y[4],y[5]);
                Blaze3DVec Fieldv(3UL, Field);
                Blaze8DVec dydxv( TEvaluateRhsGivenB(yv, Fieldv) );

                for (int i = 0; i < 6; i ++)
                { dydx[i] = dydxv[i]; }
            }

        __attribute__((always_inline)) 
            Blaze8DVec TEvaluateRhsGivenB( const Blaze3DVec y,
			                                  const Blaze3DVec B ) const
            {
                Blaze8DVec dydxv;
                
                //dot product
                G4double momentum_mag_square = (y, y);
                
                G4double inv_momentum_magnitude = vdt::fast_isqrt_general( momentum_mag_square, 4);
                G4double cof = FCof()*inv_momentum_magnitude;
                
                //scalar product
                subvector(dydxv, 0UL, 3UL) = 
                    inv_momentum_magnitude*y;
		          
                //cross product
                subvector(dydxv, 3UL, 3UL) = cof*( y % B );
                
                return dydxv;
            }

    private:
        enum { G4maximum_number_of_field_components = 24 };
        T_Field *itsField;
};

