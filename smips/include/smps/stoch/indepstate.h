#ifndef SMIPS_STOCH_INDEPSTATE_H
#define SMIPS_STOCH_INDEPSTATE_H

#include "smps/parserstate.h"

namespace smps::stoch
{
    class IndepState : public smps::ParserState
    {
    public:
        bool parse(smps::Smps &smps, std::string &line) override;

        bool maybeTransition(std::unique_ptr<smps::ParserState> &state,
                             std::string &line) override;
    };
}  // namespace smps::stoch

#endif  // SMIPS_STOCH_INDEPSTATE_H
