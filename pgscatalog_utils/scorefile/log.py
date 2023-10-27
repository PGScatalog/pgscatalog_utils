import logging
import logging.handlers


def log_format() -> str:
    return '%(asctime)s %(processName)-10s %(name)s %(levelname)-8s %(message)s'


def worker_configurer(queue):
    h = logging.handlers.QueueHandler(queue)
    root = logging.getLogger()
    root.addHandler(h)
    root.setLevel(logging.DEBUG)


def listener_configurer():
    root = logging.getLogger()
    h = logging.StreamHandler()
    h.setFormatter(logging.Formatter(log_format()))
    root.addHandler(h)


def listener_process(queue, configurer):
    configurer()
    while True:
        try:
            record = queue.get()
            if record is None:  # We send this as a sentinel to tell the listener to quit.
                break
            logger = logging.getLogger(record.name)
            logger.handle(record)  # No level or filter logic applied - just do it!
        except Exception:
            break
