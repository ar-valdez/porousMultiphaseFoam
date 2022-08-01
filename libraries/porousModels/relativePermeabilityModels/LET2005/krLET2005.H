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
    Foam::relativePermeabilityModels::krLET2005

Description
    Standard Corey1954 relative permeability model.

SourceFiles
    krLET2005.C

\*---------------------------------------------------------------------------*/

#ifndef krLET2005_H
#define krLET2005_H

#include "relativePermeabilityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relativePermeabilityModels
{

/*---------------------------------------------------------------------------*\
    Class krLET2005 Declaration
    \*---------------------------------------------------------------------------*/

class krLET2005
:
    public relativePermeabilityModel
{
    //- Private data
    volScalarField Smin_;
    volScalarField Smax_;

    dictionary krLET2005Coeffs_;
    volScalarField Lb_;
    volScalarField Eb_;
    volScalarField Tb_;
    volScalarField La_;
    volScalarField Ea_;
    volScalarField Ta_;

    //- end points
    volScalarField kramax_;
    volScalarField krbmax_;

    //Sb is the water saturation
    //nb is the water exponent
    //na is the oil exponent
        
public:

    //- Runtime type information
    TypeName("LET2005");

    // Constructors

    //- Construct from components
    krLET2005(const word& name,const dictionary& transportProperties,const volScalarField& Sb);

    //- Destructor
    ~krLET2005()
        {}

    // Member Functions

    //- Correct the relative permeabilities
    void correct()
        {
            Se_= (Sb_-Smin_)/(Smax_-Smin_);
            krb_ = krbmax_*(pow(Se_,Lb_) / (pow(Se_,Lb_) + Eb_*pow(scalar(1)-Se_,Tb_)));
            kra_ = kramax_*(pow(scalar(1)-Se_,La_) / (pow(scalar(1)-Se_,La_) + Ea_*pow(Se_,Ta_))); 
            dkrbdS_ = krbmax_*(scalar(1)/(Smax_- Smin_)) * (Lb_*pow(Se_,Lb_-scalar(1))*(pow(Se_,Lb_)-Eb_*pow(scalar(1)-Se_,Tb_)) + pow(Se_,Lb_)*(Eb_*Tb_*pow(scalar(1)-Se_,Tb_-scalar(1))-Lb_*pow(Se_,Lb_-scalar(1)))) / (pow((pow(Se_,Lb_) + Eb_*pow(scalar(1)-Se_,Tb_)),scalar(2)));
            dkradS_ = kramax_*(scalar(1)/(Smax_- Smin_))*(pow(scalar(1)-Se_,La_) * (La_*pow(scalar(1)-Se_,La_-scalar(1)) - Ea_*Ta_*pow(Se_,Ta_-scalar(1)))-La_*pow(scalar(1)-Se_,La_-scalar(1)) * (pow(scalar(1)-Se_,La_) + Ea_*pow(Se_,Ta_))) / (pow((pow(scalar(1)-Se_,La_) + Ea_*pow(Se_,Ta_)),scalar(2)));
        }
    void correcta()
        {
            Se_= (Se_-Smin_)/(Smax_-Smin_);
            kra_ = kramax_*(pow(scalar(1)-Se_,La_) / (pow(scalar(1)-Se_,La_) + Ea_*pow(Se_,Ta_))); 
            dkradS_ = kramax_*(scalar(1)/(Smax_- Smin_))*(pow(scalar(1)-Se_,La_) * (La_*pow(scalar(1)-Se_,La_-scalar(1)) - Ea_*Ta_*pow(Se_,Ta_-scalar(1)))-La_*pow(scalar(1)-Se_,La_-scalar(1)) * (pow(scalar(1)-Se_,La_) + Ea_*pow(Se_,Ta_))) / (pow((pow(scalar(1)-Se_,La_) + Ea_*pow(Se_,Ta_)),scalar(2)));
        }
    void correctb()
        {
            Se_= (Sb_-Smin_)/(Smax_-Smin_);
            krb_ = krbmax_*(pow(Se_,Lb_) / (pow(Se_,Lb_) + Eb_*pow(scalar(1)-Se_,Tb_)));
            dkrbdS_ = krbmax_*(scalar(1)/(Smax_- Smin_)) * (Lb_*pow(Se_,Lb_-scalar(1))*(pow(Se_,Lb_)-Eb_*pow(scalar(1)-Se_,Tb_)) + pow(Se_,Lb_)*(Eb_*Tb_*pow(scalar(1)-Se_,Tb_-scalar(1))-Lb_*pow(Se_,Lb_-scalar(1)))) / (pow((pow(Se_,Lb_) + Eb_*pow(scalar(1)-Se_,Tb_)),scalar(2)));
        }
    void correctb(const label& celli)
        {
            scalar Se = (Sb_[celli]-Smin_[celli])/(Smax_[celli]-Smin_[celli]);
            krb_[celli] = krbmax_[celli]*(pow(Se,Lb_[celli]) / (pow(Se,Lb_[celli]) + Eb_[celli]*pow(scalar(1)-Se,Tb_[celli])));
            dkrbdS_ = krbmax_[celli]*(scalar(1)/(Smax_[celli]- Smin_[celli])) * (Lb_[celli]*pow(Se_,Lb_[celli]-scalar(1))*(pow(Se_,Lb_[celli])-Eb_[celli]*pow(scalar(1)-Se_,Tb_[celli])) + pow(Se_,Lb_[celli])*(Eb_[celli]*Tb_[celli]*pow(scalar(1)-Se_,Tb_[celli]-scalar(1))-Lb_[celli]*pow(Se_,Lb_[celli]-scalar(1)))) / (pow((pow(Se,Lb_[celli]) + Eb_[celli]*pow(scalar(1)-Se,Tb_[celli])),scalar(2)));
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relativePermeabilityModels

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
