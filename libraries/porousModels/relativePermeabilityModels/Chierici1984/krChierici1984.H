/*---------------------------------------------------------------------------*\
  		  _______  ____    ____  ________  
 		 |_   __ \|_   \  /   _||_   __  | 
   		   | |__) | |   \/   |    | |_ \_| 
   		   |  ___/  | |\  /| |    |  _|    
    		  _| |_    _| |_\/_| |_  _| |_     
   		 |_____|  |_____||_____||_____|    
   	     Copyright (C) Toulouse INP, Pierre Horgue

License
    This file is part of porousMultiphaseFoam, an extension of OpenFOAM
    developed by Pierre Horgue (phorgue@imft.fr) and dedicated to multiphase 
    flows through porous media.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::relativePermeabilityModels::krChierici1984

Description
    Standard Corey1954 relative permeability model.

SourceFiles
    krChierici1984.C

\*---------------------------------------------------------------------------*/

#ifndef krChierici1984_H
#define krChierici1984_H

#include "relativePermeabilityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{

/*---------------------------------------------------------------------------*\
    Class krChierici1984 Declaration
    \*---------------------------------------------------------------------------*/

class krChierici1984
:
    public relativePermeabilityModel
{
    //- Private data
    volScalarField Smin_;
    volScalarField Smax_;

    dictionary krChierici1984Coeffs_;
    volScalarField cB_;
    volScalarField cM_;
    volScalarField cA_;
    volScalarField cL_;

    //- end points
    volScalarField kramax_;
    volScalarField krbmax_;

    //Sb is the water saturation
        
public:

    //- Runtime type information
    TypeName("Chierici1984");

    // Constructors

    //- Construct from components
    krChierici1984(const word& name,const dictionary& transportProperties,const volScalarField& Sb);

    //- Destructor
    ~krChierici1984()
        {}

    // Member Functions

    //- Correct the relative permeabilities
    void correct()
        {
            if(Sb_ < Smin_)
            {
                Se_=scalar(0);
            }
            else if(Sb_ > Smax_)
            {
                Se_=scalar(1);
            }
            else
            {
                Se_= (Sb_-Smin_)/(Smax_-Sb_);
            }
            kra_ = kramax_ * exp(-cA_*pow(Se_,cL_));
            krb_ = krbmax_ * exp(-cB_*pow(Se_,-cM_));
            dkradS_ = -cA_*cL_*pow(Se_,cL_-scalar(1))*kra_*     (Smax_-Smin_)/(pow(Smax_-Se_,scalar(2)));
            dkrbdS_ = cB_*cM_*pow(Se_,-cM_-scalar(1))*krb_*     (Smax_-Smin_)/(pow(Smax_-Se_,scalar(2)));
        }
    void correcta()
        {
            if(Sb_ < Smin_)
            {
                Se_=scalar(0);
            }
            else if(Sb_ > Smax_)
            {
                Se_=scalar(1);
            }
            else
            {
                Se_= (Sb_-Smin_)/(Smax_-Sb_);
            }
            kra_ = kramax_ * exp(-cA_*pow(Se_,cL_));
            dkradS_ = -cA_*cL_*pow(Se_,cL_-scalar(1))*kra_*(Smax_-Smin_)/(pow(Smax_-Se_,scalar(2)));
        }
    void correctb()
        {
            if(Sb_ < Smin_)
            {
                Se_=scalar(0);
            }
            else if(Sb_ > Smax_)
            {
                Se_=scalar(1);
            }
            else
            {
                Se_= (Sb_-Smin_)/(Smax_-Sb_);
            }
            krb_ = krbmax_ * exp(-cB_*pow(Se_,-cM_));
            dkrbdS_ = cB_*cM_*pow(Se_,-cM_-scalar(1))*krb_*(Smax_-Smin_)/(pow(Smax_-Se_,scalar(2)));
        }
    void correctb(const label& celli)
        {
            if(Sb_[celli] < Smin_[celli])
            {
                scalar Se =scalar(0);
                krb_[celli] = krbmax_[celli] * exp(-cB_[celli]*pow(Se,-cM_[celli]));
            }
            else if(Sb_[celli] > Smax_[celli])
            {
                scalar Se =scalar(1);
                krb_[celli] = krbmax_[celli] * exp(-cB_[celli]*pow(Se,-cM_[celli]));
            }
            else
            {
                scalar Se = (Sb_[celli]-Smin_[celli])/(Smax_[celli]-Sb_[celli]);
                krb_[celli] = krbmax_[celli] * exp(-cB_[celli]*pow(Se,-cM_[celli]));
            }
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativePermeabilityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
