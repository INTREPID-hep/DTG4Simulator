# CommandLineParser

Lightweight singleton for parsing command-line arguments in DTSim. Provides global access for Geant4 classes where constructor arguments cannot be modified.

## Quick Start

```cpp
// In main()
auto* parser = DTSim::CommandLineParser::Instance();

// Register options (support aliases: "-m,--macro")
parser->AddFlag("-b,--batch", "Run in batch mode");
parser->AddOption("-m,--macro", "Macro file", false, "vis.mac");
parser->AddOption("-o,--output", "Output file", true);  // required

// Parse
if (parser->Parse(argc, argv) != 0) {
    DTSim::CommandLineParser::DeleteInstance();
    return 1;
}

// Use anywhere
G4bool batch = parser->HasFlag("-b");
G4String macro = parser->GetOption("-m");

// Cleanup at end
DTSim::CommandLineParser::DeleteInstance();
```

## API

### Registration
```cpp
void AddFlag(const G4String& flag, const G4String& description);
// Flags: "-b" or "-b,--batch"

void AddOption(const G4String& option, const G4String& description,
               G4bool required, const G4String& defaultValue);
// required=true → no default allowed
// required=false → default required
```

### Parsing
```cpp
G4int Parse(G4int argc, char** argv);
// Returns: 0=success, 1=help, -1=error
```

### Access
```cpp
G4bool HasFlag(const G4String& flag) const;
G4String GetOption(const G4String& option) const;
```

## Usage in Geant4 Classes

```cpp
// DetectorConstruction.cc
auto* parser = DTSim::CommandLineParser::Instance();
if (parser->HasFlag("-B")) {
    // Enable magnetic field
}
```

## Features

- **Aliases**: `-b,--batch` syntax for short/long forms
- **Validation**: Enforces required/default logic at registration
- **Auto-help**: Built-in `-h,--help` 
- **Error detection**: Unknown options, missing values, missing required options

## Example

```bash
# Help
./exampleDTSim --help

# Usage
./exampleDTSim -b -m run.mac -E 50
./exampleDTSim --batch --macro run.mac --energy 50
```

## Notes

- Singleton allows access from any Geant4 class
- Thread-safe (accessed before MT initialization)
- Call `DeleteInstance()` at end of main
