#include "DTSimLogger.hh"

namespace DTSim {

Logger* Logger::fInstance = nullptr;

Logger* Logger::Instance() {
    if (!fInstance) fInstance = new Logger();
    return fInstance;
}

Logger::Logger() : fGlobalLevel(LogLevel::INFO), fMessenger(nullptr) {
    DefineCommands();
}

Logger::~Logger() { delete fMessenger; }

bool Logger::ShouldLog(const G4String& module, LogLevel level) const {
    LogLevel target = fGlobalLevel;
    auto it = fModuleLevels.find(module);
    if (it != fModuleLevels.end()) target = it->second;
    return level <= target;
}

void Logger::DefineCommands() {
    fMessenger = new G4GenericMessenger(this, "/DTSim/log/", "Logging control");
    
    // Global level: /DTSim/log/level
    fMessenger->DeclareMethod("level", &Logger::SetGlobalLevel)
        .SetGuidance("Set global log level (0:Error, 1:Warn, 2:Info, 3:Debug)")
        .SetParameterName("lvl", false);

    // Module specific: /DTSim/log/moduleLevel
    fMessenger->DeclareMethod("moduleLevel", &Logger::SetModuleLevel)
        .SetGuidance("Set log level for a specific module")
        .SetParameterName("module lvl", false);
}

void Logger::SetGlobalLevel(G4int level) { 
    fGlobalLevel = static_cast<LogLevel>(level); 
}

void Logger::SetModuleLevel(G4String module, G4int level) {
    fModuleLevels[module] = static_cast<LogLevel>(level);
}

LogProxy::LogProxy(const G4String& module, LogLevel level) {
    fActive = Logger::Instance()->ShouldLog(module, level);
    if (fActive) {
        G4String head = (level == LogLevel::DEBUG) ? "[#] DEBUG" : 
                        (level == LogLevel::WARN)  ? "[?] WARN "  :
                        (level == LogLevel::ERROR) ? "[!] ERROR" : "[i] INFO ";
        G4cout << head << " (" << module << "): ";
    }
}

} // namespace DTSim