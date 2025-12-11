// Simple command-line argument parser for DTSim
// Lightweight implementation without unnecessary complexity

#ifndef DTSIM_COMMANDLINEPARSER_HH
#define DTSIM_COMMANDLINEPARSER_HH

#include "globals.hh"
#include <map>
#include <string>

namespace DTSim
{

/// \brief Simple command-line argument parser (Singleton)
/// Provides command-line parsing with automatic help generation
/// 
/// Singleton allows access from Geant4 user classes (DetectorConstruction, etc.)
/// where constructor arguments cannot be modified.
/// 
/// Usage in main():
///   auto* parser = CommandLineParser::Instance();
///   parser->AddOption("-m", "Macro file", true);
///   parser->AddFlag("-B", "Include magnetic field");
///   if (parser->Parse(argc, argv)) return 1;
///   
/// Usage in DetectorConstruction:
///   auto* parser = CommandLineParser::Instance();
///   if (parser->HasFlag("-B")) { EnableMagneticField(); }
class CommandLineParser
{
  public:
    /// Get singleton instance (creates if doesn't exist)
    static CommandLineParser* Instance();
    
    /// Delete singleton instance (call at end of main)
    static void DeleteInstance();
    
    /// Prevent copying
    CommandLineParser(const CommandLineParser&) = delete;
    CommandLineParser& operator=(const CommandLineParser&) = delete;
    
    /// Add a flag option (no value required, e.g., -b, --help)
    /// \param flag The flag name (e.g., "-b", "--help")
    /// \param description Help text description
    void AddFlag(const G4String& flag, const G4String& description = "");
    
    /// Add an option that requires a value (e.g., -m file.mac)
    /// \param option The option name (e.g., "-m")
    /// \param description Help text description
    /// \param required If true, option must be provided
    /// \param defaultValue Default value if not provided (only for non-required)
    void AddOption(const G4String& option, const G4String& description = "",
                   G4bool required = false, const G4String& defaultValue = "");
    
    /// Parse command-line arguments
    /// \return 0 on success, 1 if help requested, -1 on error
    G4int Parse(G4int argc, char** argv);
    
    /// Check if a flag was provided
    G4bool HasFlag(const G4String& flag) const;
    
    /// Get the value of an option
    /// \return The option value, or empty string if not provided
    G4String GetOption(const G4String& option) const;
    
    /// Get the program name (argv[0])
    const G4String& GetProgramName() const { return fProgramName; }
    
    /// Print help message
    void PrintHelp() const;

  private:
    /// Private constructor for singleton
    CommandLineParser();
    
    /// Destructor
    ~CommandLineParser() = default;
    
    struct OptionInfo {
      G4String description;
      G4bool requiresValue;
      G4bool isRequired;
      G4String defaultValue;
      G4bool wasProvided;
      G4String value;
      G4String primaryName;  // The main option name
      G4String alias;        // Alternative name (e.g., long form)
    };
    
    static CommandLineParser* fInstance;
    std::map<G4String, OptionInfo> fOptions;
    G4String fProgramName;
    
    /// Check if string looks like an option flag
    G4bool IsFlag(const char* str) const;
};

}  // namespace DTSim

#endif  // DTSIM_COMMANDLINEPARSER_HH
