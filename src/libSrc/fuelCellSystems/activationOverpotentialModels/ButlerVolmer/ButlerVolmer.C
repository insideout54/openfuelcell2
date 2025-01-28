/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
\*---------------------------------------------------------------------------*/

#include "ButlerVolmer.H"
#include "phaseModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::activationOverpotentialModels::ButlerVolmer<Thermo>::ButlerVolmer
(
    const phaseModel& phase,
    const dictionary& dict
)
:
    ActivationOverpotentialModel<Thermo>(phase, dict)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::activationOverpotentialModels::ButlerVolmer<Thermo>::~ButlerVolmer()
{}

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //
template<class Thermo>
void Foam::activationOverpotentialModels::ButlerVolmer<Thermo>::correct()
{
    // Catalyst properties (hardcoded)
    const scalar mPt = 0.5;    // Pt loading [mg/cm^2], convert if necessary
    const scalar As = 60.0;    // Specific surface area [m^2/mg]

    // Compute specific active area Av
    const scalar Av = mPt * As;

    // Create Nernst model
    this->createNernst();
    this->nernst_->correct();

    // Access relevant regions
    const regionType& fluidPhase = this->region
    (
        word(this->regions_.subDict("fluid").lookup("name"))
    );
    const regionType& electronPhase = this->region
    (
        word(this->regions_.subDict("electron").lookup("name"))
    );
    const regionType& ionPhase = this->region
    (
        word(this->regions_.subDict("ion").lookup("name"))
    );

    // Source/sink fields
    scalarField& SE = const_cast<volScalarField&>
    (
        electronPhase.template lookupObject<volScalarField>("J")
    );
    scalarField& SI = const_cast<volScalarField&>
    (
        ionPhase.template lookupObject<volScalarField>("J")
    );

    // Potential fields
    const scalarField& phiE = electronPhase.template
        lookupObject<volScalarField>(word(this->phiNames_["electron"]));
    const scalarField& phiI = ionPhase.template
        lookupObject<volScalarField>(word(this->phiNames_["ion"]));

    // Overpotential fields
    scalarField& eta = this->eta_;
    scalarField& nernst = this->nernst_()();
    scalarField& j = this->j_;

    // Local T, species, etc.
    const scalarField& T = this->thermo_.T();
    const scalarField& s = this->phase_;

    // Mixture mole fraction and density
    const scalarField W(this->thermo_.W() / 1000.0);
    const scalarField& rho = this->phase_.thermo().rho();

    scalarField coeff(this->phase_.mesh().nCells(), 1.0);
    forAllConstIter(dictionary, this->species_, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& d0 = iter().dict();

            scalar ksi = readScalar(d0.lookup("ksi"));
            scalar cRef = readScalar(d0.lookup("cRef"));

            const scalarField& X = this->phase_.X(name);

            coeff *= Foam::pow(X * rho / W / cRef, ksi);
        }
    }

    // Catalyst zone information
    label znId = fluidPhase.cellZones().findZoneID(this->zoneName_);
    const labelList& cells = fluidPhase.cellZones()[znId];

    // Electron transfer sign
    scalar ePerMol = this->nernst_->rxnList()["e"];
    scalar sign = ePerMol / mag(ePerMol);

    // Accumulate total current
    scalar Rj(0.0);

    forAll(cells, cellI)
    {
        label fluidId = cells[cellI];
        label electronId = electronPhase.cellMap()[fluidPhase.cellMapIO()[fluidId]];
        label ionId = ionPhase.cellMap()[fluidPhase.cellMapIO()[fluidId]];

        // Overpotential relaxation
        eta[fluidId] = 
        (
            eta[fluidId] * (scalar(1) - this->relax_)
          + (phiE[electronId] - phiI[ionId] - nernst[fluidId]) * this->relax_
        );

        // Clamping the exponents for stability
        scalar localEtaAn = ePerMol * this->alpha_ * F * eta[fluidId] / (Rgas * T[fluidId]);
        scalar localEtaCat = -ePerMol * (scalar(1) - this->alpha_) * F * eta[fluidId] / (Rgas * T[fluidId]);

        localEtaAn = Foam::clamp(localEtaAn, -25.0, 25.0);
        localEtaCat = Foam::clamp(localEtaCat, -25.0, 25.0);

        // Reactant concentration clamp
        scalar sLocal = Foam::max(s[fluidId], VSMALL);

        // Compute current density using Butler–Volmer with Av
        scalar expAn = Foam::exp(localEtaAn);
        scalar expCat = Foam::exp(localEtaCat);

        scalar jLocal = Av * this->j0_.value() * coeff[fluidId] * Foam::pow(sLocal, this->gamma_)
                       * (expAn - expCat);

        j[fluidId] = Foam::max(jLocal, scalar(0));

        // Source/sink for electron & ion phases
        SE[electronId] = -sign * j[fluidId];
        SI[ionId] = -SE[electronId];

        Rj += fluidPhase.V()[fluidId] * j[fluidId];
    }

    reduce(Rj, sumOp<scalar>());

    Info << "Total current (A) at " << this->zoneName_ << ": " << Rj << endl;
}

// ************************************************************************* //
