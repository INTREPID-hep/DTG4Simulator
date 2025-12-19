#ifndef DTSimLogger_hh
#define DTSimLogger_hh 1

#include "G4ios.hh"
#include "G4String.hh"
#include "G4GenericMessenger.hh"
#include <map>

namespace DTSim {

enum class LogLevel { ERROR = 0, WARN = 1, INFO = 2, DEBUG = 3 };

class Logger {
public:
    static Logger* Instance();
    
    // Verbosity Control
    void SetGlobalLevel(G4int level);
    void SetModuleLevel(G4String module, G4int level);
    bool ShouldLog(const G4String& module, LogLevel level) const;

private:
    Logger();
    ~Logger();
    void DefineCommands();

    static Logger* fInstance;
    LogLevel fGlobalLevel;
    std::map<G4String, LogLevel> fModuleLevels;
    G4GenericMessenger* fMessenger;
};

class LogProxy {
public:
    LogProxy(const G4String& module, LogLevel level);
    ~LogProxy() = default;

    template<typename T>
    LogProxy& operator<<(const T& msg) {
        if (fActive) G4cout << msg;
        return *this;
    }

    LogProxy& operator<<(std::ostream& (*f)(std::ostream&)) {
        if (fActive) G4cout << f;
        return *this;
    }

private:
    bool fActive;
};

// Simplified usage macros
#define LogError(mod) DTSim::LogProxy(mod, DTSim::LogLevel::ERROR)
#define LogWarn(mod)  DTSim::LogProxy(mod, DTSim::LogLevel::WARN)
#define LogInfo(mod)  DTSim::LogProxy(mod, DTSim::LogLevel::INFO)
#define LogDebug(mod) DTSim::LogProxy(mod, DTSim::LogLevel::DEBUG)

} // namespace DTSim
#endif