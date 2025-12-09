#include "DTSimUtils.hh"
#include "G4String.hh"
#include <string>

namespace DTSim
{

G4int ExtractIntAfterToken(const G4String& str, const G4String& token, size_t startPos)
{
    // Find the token in the string
    size_t pos = str.find(token, startPos);
    if (pos == std::string::npos) {
        return -999;  // Token not found
    }
    
    // Move past the token
    pos += token.length();
    
    // Find the end of the number (next underscore or end of string)
    size_t endPos = str.find("_", pos);
    
    // Extract the substring containing the number
    G4String numberStr = str.substr(pos, endPos - pos);
    
    // Convert to integer (handles negative numbers)
    try {
        return std::stoi(numberStr);
    } catch (const std::exception& e) {
        return -999;  // Parsing failed
    }
}

}
