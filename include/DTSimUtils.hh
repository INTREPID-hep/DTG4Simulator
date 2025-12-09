#ifndef DTSimUtils_hh
#define DTSimUtils_hh 1

#include "G4Types.hh"
#include "G4String.hh"

namespace DTSim
{

/**
 * @brief Extract an integer value that appears after a specific token in a string
 * 
 * @param str The string to search in (e.g., volume name)
 * @param token The token to search for (e.g., "_W", "_Sec", "_St")
 * @param startPos Optional starting position for the search
 * @return The extracted integer value, or -999 if token not found or parsing fails
 * 
 * Example: ExtractIntAfterToken("Station_W2_Sec5_St3", "_W") returns 2
 *          ExtractIntAfterToken("Station_W-1_Sec5_St3", "_W") returns -1
 */
G4int ExtractIntAfterToken(const G4String& str, const G4String& token, 
                          size_t startPos = 0);

}

#endif
