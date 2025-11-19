// Simple command-line argument parser for DTSim

#include "CommandLineParser.hh"
#include "G4Exception.hh"
#include <iomanip>
#include <cstring>

using namespace DTSim;

CommandLineParser* CommandLineParser::fInstance = nullptr;

namespace {
  /// Helper function to split comma-separated option names
  /// E.g., "-b,--batch" → primary="-b", alias="--batch"
  void SplitOptionNames(const G4String& input, G4String& primary, G4String& alias)
  {
    size_t pos = input.find(',');
    if (pos == std::string::npos) {
      primary = input;
      alias = "";
      return;
    }
    
    primary = input.substr(0, pos);
    alias = input.substr(pos + 1);
    
    // Trim whitespace from alias
    while (!alias.empty() && alias.front() == ' ') alias.erase(0, 1);
    while (!alias.empty() && alias.back() == ' ') alias.pop_back();
  }
}


CommandLineParser* CommandLineParser::Instance()
{
  if (!fInstance) {
    fInstance = new CommandLineParser();
  }
  return fInstance;
}


void CommandLineParser::DeleteInstance()
{
  if (fInstance) {
    delete fInstance;
    fInstance = nullptr;
  }
}


CommandLineParser::CommandLineParser()
  : fProgramName("")
{
  // Register built-in help flags with alias
  AddFlag("-h,--help", "Print this help message");
}


void CommandLineParser::AddFlag(const G4String& flag, const G4String& description)
{
  // Parse comma-separated names: "-b,--batch" or just "-b"
  G4String primaryName, aliasName;
  SplitOptionNames(flag, primaryName, aliasName);
  
  OptionInfo info;
  info.description = description;
  info.requiresValue = false;
  info.isRequired = false;
  info.wasProvided = false;
  info.primaryName = primaryName;
  info.alias = aliasName;
  
  // Register primary name
  fOptions[primaryName] = info;
  
  // Register alias if exists
  if (!info.alias.empty()) {
    fOptions[info.alias] = info;
  }
}


void CommandLineParser::AddOption(const G4String& option, const G4String& description,
                                  G4bool required, const G4String& defaultValue)
{
  // Validate: required options shouldn't have defaults
  if (required && !defaultValue.empty()) {
    G4ExceptionDescription desc;
    desc << "Option '" << option << "' is marked required but has default value '" 
         << defaultValue << "'. Required options cannot have defaults.";
    G4Exception("CommandLineParser::AddOption", "CLI003", FatalException, desc);
  }
  
  // Validate: optional options should have defaults
  if (!required && defaultValue.empty()) {
    G4ExceptionDescription desc;
    desc << "Option '" << option << "' is optional but has no default value. "
         << "Either make it required (true) or provide a default value.";
    G4Exception("CommandLineParser::AddOption", "CLI004", FatalException, desc);
  }
  
  // Parse comma-separated names: "-m,--macro" or just "-m"
  G4String primaryName, aliasName;
  SplitOptionNames(option, primaryName, aliasName);
  
  OptionInfo info;
  info.description = description;
  info.requiresValue = true;
  info.isRequired = required;
  info.defaultValue = defaultValue;
  info.wasProvided = false;
  info.value = defaultValue;
  info.primaryName = primaryName;
  info.alias = aliasName;
  
  // Register primary name
  fOptions[primaryName] = info;
  
  // Register alias if exists
  if (!info.alias.empty()) {
    fOptions[info.alias] = info;
  }
}


G4bool CommandLineParser::IsFlag(const char* str) const
{
  return str != nullptr && str[0] == '-';
}


G4int CommandLineParser::Parse(G4int argc, char** argv)
{
  if (argc < 1) return -1;
  
  fProgramName = argv[0];
  
  // Parse arguments
  for (G4int i = 1; i < argc; ++i) {
    G4String arg = argv[i];
    
    // Check if this is a registered option
    auto it = fOptions.find(arg);
    if (it == fOptions.end()) {
      G4cerr << "Error: Unknown option '" << arg << "'" << G4endl;
      PrintHelp();
      return -1;
    }
    
    OptionInfo& opt = it->second;
    opt.wasProvided = true;
    
    // If option requires a value, get the next argument
    if (opt.requiresValue) {
      if (i + 1 >= argc) {
        G4cerr << "Error: Option '" << arg << "' requires a value" << G4endl;
        return -1;
      }
      
      ++i;
      if (IsFlag(argv[i])) {
        G4cerr << "Error: Option '" << arg << "' requires a value, got flag '" 
               << argv[i] << "'" << G4endl;
        return -1;
      }
      
      opt.value = argv[i];
    }
  }
  
  // Check if help was requested
  if (HasFlag("-h") || HasFlag("--help")) {
    PrintHelp();
    return 1;
  }
  
  // Check required options
  for (const auto& pair : fOptions) {
    const OptionInfo& opt = pair.second;
    if (opt.isRequired && !opt.wasProvided) {
      G4cerr << "Error: Required option '" << pair.first << "' not provided" << G4endl;
      PrintHelp();
      return -1;
    }
  }
  
  return 0;
}


G4bool CommandLineParser::HasFlag(const G4String& flag) const
{
  auto it = fOptions.find(flag);
  if (it == fOptions.end()) {
    G4ExceptionDescription desc;
    desc << "Flag '" << flag << "' was not registered. Call AddFlag() first.";
    G4Exception("CommandLineParser::HasFlag", "CLI001", FatalException, desc);
    return false;
  }
  
  return it->second.wasProvided;
}


G4String CommandLineParser::GetOption(const G4String& option) const
{
  auto it = fOptions.find(option);
  if (it == fOptions.end()) {
    G4ExceptionDescription desc;
    desc << "Option '" << option << "' was not registered. Call AddOption() first.";
    G4Exception("CommandLineParser::GetOption", "CLI002", FatalException, desc);
    return "";
  }
  
  return it->second.value;
}


void CommandLineParser::PrintHelp() const
{
  G4cout << "\nUsage: " << fProgramName << " [OPTIONS]\n" << G4endl;
  G4cout << "Options:" << G4endl;
  
  // Calculate max width for alignment (including aliases)
  size_t maxWidth = 0;
  for (const auto& pair : fOptions) {
    const OptionInfo& opt = pair.second;
    // Only calculate width for primary names
    if (pair.first != opt.primaryName) continue;
    
    size_t width = opt.primaryName.length();
    if (!opt.alias.empty()) {
      width += 2 + opt.alias.length();  // ", --alias"
    }
    if (opt.requiresValue) {
      width += 8;  // " <value>"
    }
    if (width > maxWidth) maxWidth = width;
  }
  
  // Print each option (only primaries, aliases shown inline)
  for (const auto& pair : fOptions) {
    const G4String& name = pair.first;
    const OptionInfo& opt = pair.second;
    
    // Skip if this is an alias (not the primary name)
    if (name != opt.primaryName) continue;
    
    // Build option display string
    G4String line = "  " + opt.primaryName;
    
    // Add alias if exists
    if (!opt.alias.empty()) {
      line += ", " + opt.alias;
    }
    
    // Add value indicator for options
    if (opt.requiresValue) {
      line += " <value>";
    }
    
    // Padding for alignment
    while (line.length() < maxWidth + 4) {
      line += " ";
    }
    
    // Add description
    line += opt.description;
    
    // Add default value or required indicator
    if (opt.requiresValue) {
      if (opt.isRequired) {
        line += " [required]";
      } else if (!opt.defaultValue.empty()) {
        line += " [default: " + opt.defaultValue + "]";
      }
    }
    
    G4cout << line << G4endl;
  }
  
  G4cout << G4endl;
}
