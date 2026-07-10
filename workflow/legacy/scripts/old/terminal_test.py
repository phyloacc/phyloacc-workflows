import sys
import time
import threading
import queue
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from multiprocessing import Manager

def log(msg):
    sys.stdout.write(msg + "\n")
    sys.stdout.flush()

def test_carriage_return():
    log("TEST 1: carriage-return single-line update")
    for i in range(21):
        sys.stdout.write(f"\rprogress {i:02d}/20")
        sys.stdout.flush()
        time.sleep(0.1)
    sys.stdout.write("\n")
    sys.stdout.flush()

def test_ansi_rewrite():
    log("TEST 2: ANSI multi-line rewrite")
    sys.stdout.write("line A\n")
    sys.stdout.write("line B\n")
    sys.stdout.flush()
    time.sleep(0.5)

    for i in range(10):
        sys.stdout.write("\x1b[2F")   # move up 2 lines
        sys.stdout.write("\x1b[2K")
        sys.stdout.write(f"line A updated {i}\n")
        sys.stdout.write("\x1b[2K")
        sys.stdout.write(f"line B updated {i}\n")
        sys.stdout.flush()
        time.sleep(0.2)

def test_threaded_rewrite():
    log("TEST 3: background-thread rewrite")
    stop = threading.Event()

    def painter():
        tick = 0
        while not stop.is_set():
            sys.stdout.write("\x1b[2K")
            sys.stdout.write(f"\rthread progress tick={tick}")
            sys.stdout.flush()
            tick += 1
            time.sleep(0.2)

    t = threading.Thread(target=painter, daemon=True)
    t.start()

    for i in range(10):
        time.sleep(0.3)

    stop.set()
    t.join()
    sys.stdout.write("\n")
    sys.stdout.flush()

def render_batch_progress(progress_queue, stop_event):
    previous_line_count = 0
    batch_state = {}
    done_after_render = set()

    while not stop_event.is_set() or batch_state:
        try:
            event_type, batch_num, state = progress_queue.get(timeout=0.2)
        except queue.Empty:
            continue

        if event_type in ("start", "update"):
            batch_state[batch_num] = state
        elif event_type == "done":
            batch_state[batch_num] = state
            done_after_render.add(batch_num)

        try:
            while True:
                event_type, batch_num, state = progress_queue.get_nowait()
                if event_type in ("start", "update"):
                    batch_state[batch_num] = state
                    done_after_render.discard(batch_num)
                elif event_type == "done":
                    batch_state[batch_num] = state
                    done_after_render.add(batch_num)
        except queue.Empty:
            pass

        snapshot = sorted(batch_state.items())
        if snapshot:
            if previous_line_count:
                sys.stdout.write(f"\x1b[{previous_line_count}F")
            for batch_num, (ok_count, warn_count, _processed, total_count) in snapshot:
                sys.stdout.write("\x1b[2K")
                sys.stdout.write(f"[ INFO ] Batch {batch_num}: {ok_count} successful / {warn_count} warn / {total_count} total\n")
            sys.stdout.flush()
            previous_line_count = len(snapshot)
            for batch_num in list(done_after_render):
                batch_state.pop(batch_num, None)
                done_after_render.discard(batch_num)

    if previous_line_count:
        sys.stdout.write(f"\x1b[{previous_line_count}F")
        for _ in range(previous_line_count):
            sys.stdout.write("\x1b[2K\n")
        sys.stdout.write(f"\x1b[{previous_line_count}F")
        sys.stdout.flush()

def forward_progress_events(manager_queue, render_queue, debug_queue, stop_event):
    while not stop_event.is_set():
        try:
            event = manager_queue.get(timeout=0.2)
        except queue.Empty:
            continue
        render_queue.put(event)
        debug_queue.put(event)

def worker_batch(batch_num, total_count, manager_queue):
    manager_queue.put(("start", batch_num, (0, 0, 0, total_count)))
    ok_count = 0
    warn_count = 0
    for i in range(1, total_count + 1):
        if i % 5 == 0:
            warn_count += 1
        else:
            ok_count += 1
        if i % 10 == 0 or i == total_count:
            manager_queue.put(("update", batch_num, (ok_count, warn_count, i, total_count)))
        time.sleep(0.03)
    manager_queue.put(("done", batch_num, (ok_count, warn_count, total_count, total_count)))
    return batch_num

def test_parent_owned_progress():
    log("TEST 4: parent-owned queue progress with worker processes")
    with Manager() as manager:
        manager_queue = manager.Queue()
        render_queue = queue.Queue()
        debug_queue = queue.Queue()
        render_stop = threading.Event()
        forward_stop = threading.Event()

        forwarder = threading.Thread(
            target=forward_progress_events,
            args=(manager_queue, render_queue, debug_queue, forward_stop),
            daemon=True,
        )
        renderer = threading.Thread(
            target=render_batch_progress,
            args=(render_queue, render_stop),
            daemon=True,
        )
        forwarder.start()
        renderer.start()

        with ProcessPoolExecutor(max_workers=3) as executor:
            futures = {executor.submit(worker_batch, batch_num, 50, manager_queue): batch_num for batch_num in range(1, 4)}
            while futures:
                done, _ = wait(list(futures.keys()), return_when=FIRST_COMPLETED)
                for future in done:
                    futures.pop(future)
                    log(f"batch {future.result()} completed")

        forward_stop.set()
        forwarder.join()
        render_stop.set()
        renderer.join()

        seen_debug = False
        try:
            while True:
                event_type, batch_num, state = debug_queue.get_nowait()
                seen_debug = True
                if event_type == "update" and state[2] == 10:
                    log(f"debug event batch {batch_num}: {state}")
        except queue.Empty:
            pass

        log(f"debug events observed = {seen_debug}")

if __name__ == "__main__":
    log(f"stdout.isatty() = {sys.stdout.isatty()}")
    log(f"stderr.isatty() = {sys.stderr.isatty()}")
    test_carriage_return()
    test_ansi_rewrite()
    test_threaded_rewrite()
    test_parent_owned_progress()
    log("done")
