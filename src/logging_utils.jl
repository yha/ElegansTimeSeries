using Dates
using LoggingExtras: EarlyFilteredLogger, FormatLogger

function progress_throttle_logger(period::Period, logger)
    history = Dict{Any, DateTime}()
    EarlyFilteredLogger(logger) do log
        if log.group == :ProgressLogging
            if !haskey(history, log.id) || (period < now() - history[log.id])
                history[log.id] = now()
                return true
            end
            return false
        end
        return true
    end
end
function logger_with_timestamps(io)
    FormatLogger(io) do io, args
        now_str = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
        println(io, "[", now_str, "] ", args.level, ": ", args.message)
        for (k,v) in args.kwargs
            println(io, " " ^ (length(now_str)+3), k, " = ", repr(v))
        end
        println(io)
    end
end
