/*---------------------------------------------------------------------------*\
  TimeAverager.H — restart-safe time-averaging helper for OpenFOAM Foundation
  Reads 'averagingStartTime' from controlDict.
\*---------------------------------------------------------------------------*/

#ifndef TimeAverager_H
#define TimeAverager_H

#include "Time.H"
#include "dictionary.H"
#include "error.H"
#include "IOstreams.H"

namespace Foam
{

class TimeAverager
{
    scalar startTime_{0.0};
    scalar currentTime_{0.0};
    scalar totalAveragingTime_{0.0};
    static constexpr scalar SMALL_TIME_ = 1e-14;

public:
    explicit TimeAverager(scalar startTime = 0.0) : startTime_(startTime) {}

    inline void setCurrentTime(const scalar t)
    {
        currentTime_ = t;
        totalAveragingTime_ = Foam::max(0.0, currentTime_ - startTime_);
    }

    inline scalar startTime() const          { return startTime_; }
    inline scalar totalAveragingTime() const { return totalAveragingTime_; }
    inline bool   active() const             { return currentTime_ >= startTime_; }

    template<class GeoField>
    inline void update(GeoField& meanField, const GeoField& instField, const scalar dt) const
    {
        if (!active())
        {
            meanField = instField; // keep in sync before averaging window
            return;
        }

        if (totalAveragingTime_ > SMALL_TIME_)
        {
            const scalar ow = (totalAveragingTime_ - dt)/totalAveragingTime_;
            const scalar nw = dt/totalAveragingTime_;
            meanField = ow*meanField + nw*instField;
			
        }
        else
        {
            meanField = instField; // first step at startTime
        }
    }
};

inline scalar readAveragingStartTime(const Time& runTime)
{
    const dictionary& controlDict = runTime.controlDict();
    const scalar t0 = controlDict.lookupOrDefault<scalar>("averagingStartTime", 0.0);
    if (t0 < 0.0)
    {
        FatalErrorInFunction
            << "averagingStartTime cannot be negative: " << t0
            << exit(FatalError);
    }
    return t0;
}

} // namespace Foam

#endif
// ************************************************************************* //

