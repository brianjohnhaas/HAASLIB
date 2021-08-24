#!/usr/bin/env python3

import time
import multiprocessing
import random
import logging


logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s : %(levelname)s : %(message)s',
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)



class MultiProcessManager:

    def __init__(self, num_parallel_processes):
        
        self.num_parallel_processes = num_parallel_processes
        self.num_running = 0
        self.num_successes = 0
        self.num_errors = 0
        self.process_list = list()

    def launch_process(self, process):

        if self.num_running >= self.num_parallel_processes:
            self.wait_for_open_slot()

        process.start()
        self.process_list.append(process)
        self.num_running += 1


    def wait_for_open_slot(self):

        while self.num_running >= self.num_parallel_processes:
            self._screen_running_processes()
            time.sleep(0.1)

    def _screen_running_processes(self):

        completed_processes = list()
            
        for i, process in enumerate(self.process_list):
            if not process.is_alive():
                completed_processes.append(i)

        if completed_processes:
            completed_processes.reverse()
            for completed_process_idx in completed_processes:
                self.num_running -= 1
                process = self.process_list[completed_process_idx]
                process.join()
                if process.exitcode == 0:
                    self.num_successes += 1
                else:
                    self.num_errors += 1
                    logger.debug("-captured a failed process")
                    
                del self.process_list[completed_process_idx]

            
                
    def wait_for_remaining_processes(self):

        logger.debug("-waiting for remaining processes")
        
        while self.num_running > 0:
            logger.debug("-waiting on {} processes".format(self.num_running))
            self._screen_running_processes()
            time.sleep(0.1)

        logger.debug("-done waiting. All processes are completed")


    def summarize_status(self):
        return("{} jobs succeeded & {} jobs failed".format(self.num_successes, self.num_errors))

        
                    
def test_mpm(num_parallel_processes=8, num_total_processes=100):

    def runner(id):
        print("running id:{}".format(id))
        x = id / (id % 10)  # should error as div-by-zero on occasion
        time.sleep(random.randint(0,10))

    mpm = MultiProcessManager(num_parallel_processes)
        
    for i in range(num_total_processes):
    
        p = multiprocessing.Process(target=runner, args=(i,))

        mpm.launch_process(p)


    mpm.wait_for_remaining_processes()
    
    



if __name__=='__main__':

    # run test
    logger.setLevel(logging.DEBUG)
    test_mpm()
