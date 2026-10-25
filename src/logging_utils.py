"""Structured logging: JSON lines to a run log file, readable lines on the console."""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from pathlib import Path

LOGGER_NAME = "extremonet"


class JsonFormatter(logging.Formatter):
    """One JSON object per line. Structured fields are passed as `extra={"ctx": {...}}`."""

    def format(self, record: logging.LogRecord) -> str:
        entry = {
            "ts": datetime.fromtimestamp(record.created, tz=timezone.utc).isoformat(timespec="milliseconds"),
            "level": record.levelname,
            "logger": record.name,
            "event": record.getMessage(),
        }
        ctx = getattr(record, "ctx", None)
        if ctx:
            entry["ctx"] = ctx
        if record.exc_info:
            entry["exception"] = self.formatException(record.exc_info)
        return json.dumps(entry, default=str)


class ConsoleFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        line = f"{record.levelname:<7} {record.getMessage()}"
        ctx = getattr(record, "ctx", None)
        if ctx:
            line += "  " + " ".join(f"{k}={v}" for k, v in ctx.items())
        return line


def get_logger(name: str | None = None) -> logging.Logger:
    return logging.getLogger(LOGGER_NAME if not name else f"{LOGGER_NAME}.{name}")


def configure_logging(level: str = "INFO", log_file: Path | None = None) -> logging.Logger:
    """Reset the package logger's handlers; safe to call once per pipeline run."""
    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(level)
    logger.propagate = False
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()

    console = logging.StreamHandler()
    console.setFormatter(ConsoleFormatter())
    logger.addHandler(console)

    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setFormatter(JsonFormatter())
        logger.addHandler(file_handler)
    return logger
