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
    Foam::relativePermeabilityModels::krCorey1954

Description
    Standard Corey1954 relative permeability model.

SourceFiles
    krCorey1954.C

\*---------------------------------------------------------------------------*/

#ifndef krCorey1954_H
#define krCorey1954_H

#include "relativePermeabilityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{

/*---------------------------------------------------------------------------*\
    Class krCorey1954 Declaration
    \*---------------------------------------------------------------------------*/

class krCorey1954
:
    public relativePermeabilityModel
{
    //- Private data
    volScalarField Smin_;
    volScalarField Smax_;

    dictionary krCorey1954Coeffs_;
    volScalarField nb_;
    volScalarField na_;

    //- end points
    volScalarField kramax_;
    volScalarField krbmax_;

    //Sb is the water saturation
    //nb is the water exponent
    //na is the oil exponent
        
public:

    //- Runtime type information
    TypeName("Corey1954");

    // Constructors

    //- Construct from components
    krCorey1954(const word& name,const dictionary& transportProperties,const volScalarField& Sb);

    //- Destructor
    ~krCorey1954()
        {}

    // Member Functions

    //- Correct the relative permeabilities
    void correct()
        {
            Se_= (Sb_-Smin_)/(Smax_-Smin_);
            krb_ = krbmax_ * pow(Se_,nb_);
            kra_ = kramax_ * pow((scalar(1)-Se_),na_);
            dkrbdS_ = krbmax_*nb_*pow(Se_,nb_-1)/(Smax_- Smin_);
            dkradS_ = -kramax_*na_*pow((scalar(1)-Se_),na_-1)/(Smax_- Smin_);
        }
    void correcta()
        {
            Se_= (Sb_-Smin_)/(Smax_-Smin_);
            kra_ = kramax_ * pow((scalar(1)-Se_),na_);
            dkradS_ = -kramax_*na_*pow((scalar(1)-Se_),na_-1)/(Smax_- Smin_);
        }
    void correctb()
        {
            Se_= (Sb_-Smin_)/(Smax_-Smin_);
            krb_ = krbmax_ * pow(Se_,nb_);
            dkrbdS_ = krbmax_*nb_*pow(Se_,nb_-1)/(Smax_- Smin_);
        }
    void correctb(const label& celli)
        {
            scalar Se = (Sb_[celli]-Smin_[celli])/(Smax_[celli]-Smin_[celli]);
            krb_[celli] = krbmax_[celli] * pow(Se,nb_[celli]);
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativePermeabilityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
