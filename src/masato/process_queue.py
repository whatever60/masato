import asyncio
from collections import deque
from threading import Event, Thread
import subprocess


class ProcessQueue:
    def __init__(self, capacity):
        self.capacity = capacity
        self.queue = deque()
        self.active = []
        self.all_done_event = Event()

    def add_process(
        self,
        args,
        stdin=None,
        stdout=None,
        stderr=None,
        shell=False,
        cwd=None,
        env=None,
    ) -> None:
        self.add_callable(
            subprocess.Popen,
            args,
            stdin=stdin,
            stdout=stdout,
            stderr=stderr,
            shell=shell,
            cwd=cwd,
            env=env,
        )

    def add_callable(self, callable_, *args, **kwargs) -> None:
        """Add a callable that starts a process."""
        if not self.queue and len(self.active) == 0:
            self.all_done_event.clear()
        self.queue.append([callable_, args, kwargs])
        self._try_to_start_process()

    def _try_to_start_process(self):
        while self.queue and len(self.active) < self.capacity:
            callable_, args, kwargs = self.queue.popleft()
            process = callable_(*args, **kwargs)
            self.active.append(process)
            thread = Thread(target=self._monitor_process, args=(process,))
            thread.start()

    def _monitor_process(self, process):
        process.wait()
        self.active.remove(process)
        if not self.active and not self.queue:
            # Set the event if there are no active processes and nothing in the queue
            self.all_done_event.set()
        else:
            self._try_to_start_process()

    def wait_for_completion(self):
        """Wait until all tasks have completed."""
        self.all_done_event.wait()


class AsyncProcessQueue:
    def __init__(self, capacity):
        self.capacity = capacity
        self.queue = deque()
        self.active = set()  # Using a set to track active tasks
        self.all_done_event = asyncio.Event()
        self.all_done_event.set()  # Initialize as set indicating no tasks are pending

    async def add_process(
        self,
        args,
        stdin=None,
        stdout=None,
        stderr=None,
        shell=False,
        cwd=None,
        env=None,
    ) -> None:
        # Wrap process creation in a callable for consistent handling
        process_callable = lambda: (
            asyncio.create_subprocess_exec(
                *args, stdin=stdin, stdout=stdout, stderr=stderr, cwd=cwd, env=env
            )
            if not shell
            else asyncio.create_subprocess_shell(
                args, stdin=stdin, stdout=stdout, stderr=stderr, cwd=cwd, env=env
            )
        )
        await self.add_callable(process_callable)

    async def add_callable(self, callable_):
        """Add a callable that starts an async process."""
        if not self.queue and not self.active:
            self.all_done_event.clear()
        self.queue.append(callable_)
        await self._try_to_start_process()

    async def _try_to_start_process(self):
        while self.queue and len(self.active) < self.capacity:
            callable_ = self.queue.popleft()
            task = asyncio.create_task(callable_())
            self.active.add(task)
            task.add_done_callback(self._task_done)

    def _task_done(self, task):
        self.active.remove(task)
        if not self.active and not self.queue:
            self.all_done_event.set()  # No more active tasks
        else:
            asyncio.create_task(
                self._try_to_start_process()
            )  # Try to start more tasks if capacity allows

    async def wait_for_completion(self):
        """Wait until all tasks have completed."""
        await self.all_done_event.wait()


async def main():
    queue = AsyncProcessQueue(capacity=2)
    commands = [
        'echo "Task 1 started"; sleep 1; echo "Task 1 completed"',
        'echo "Task 2 started"; sleep 2; echo "Task 2 completed"',
        'echo "Task 3 started"; sleep 3; echo "Task 3 completed"',
        'echo "Task 4 started"; sleep 2; echo "Task 4 completed"',
        'echo "Task 5 started"; sleep 1; echo "Task 5 completed"',
    ]
    for cmd in commands:
        await queue.add_process(cmd, shell=True)
    await queue.wait_for_completion()
    print("All processes have completed.")


if __name__ == "__main__":
    asyncio.run(main())
    # queue = ProcessQueue(capacity=2)
    # commands = [
    #     'echo "Task 1 started"; sleep 1; echo "Task 1 completed"',
    #     'echo "Task 2 started"; sleep 2; echo "Task 2 completed"',
    #     'echo "Task 3 started"; sleep 3; echo "Task 3 completed"',
    #     'echo "Task 4 started"; sleep 2; echo "Task 4 completed"',
    #     'echo "Task 5 started"; sleep 1; echo "Task 5 completed"',
    # ]
    # for cmd in commands:
    #     queue.add_process(cmd, shell=True)
    # queue.wait_for_completion()
    # print("All processes have completed.")
