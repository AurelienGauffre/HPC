Exercises/attempts to parallelize classical stuffs and algorithms using OpenMP.

## Cheatsheet for OpenMP :

➤ int omp_get_num_threads(), returns the number of threads inside a
parallel region
➤ int omp_get_thread_num(), get the thread rank inside a parallel region
➤ void omp_set_num_threads(int), sets the number of threads to be used
➤ void omp_set_nested(int), turns nested parallelism on/of

➤ #pragma omp single { ... } : The code is executed by one thread only, no guarantee which thread
Can introduce an implicit barrier at the end
➤ #pragma omp master { ... } : Code executed by the master thread, guaranteed and no implicit barrier at the
end.


➤ #pragma omp barrier :  wait for all the thread to have finished to continue, synchronization, must be encountered by all threads in a team (or none)
➤ #pragma omp ordered { a block of codes }
is another form of synchronization (in sequential order). The form
➤ #pragma omp critical { a block of codes } : only one thread at the time will work on next line ({} can be omiited for only one line)
➤ #pragma omp atomic { single assignment statement } :  same, but more efficient, only support basic operation


➤ #pragma omp for : instructs the compiler to distribute loop iterations within the team of threads (must be used inside a parallel region)


This kind of constuction can avoird using if statement on tread ID :
➤ #pragma omp parallel
{
  #pragma omp sections
  {
  #pragma omp section
  funcA ();
  #pragma omp section
  funcB ();
  #pragma omp section
  funcC ();
  }
}

waring at the end of sections, there's a call to an implicit barrier not as tasks. Section create a new parallel region. Tasks are thus more adapated for recursive implementation.
