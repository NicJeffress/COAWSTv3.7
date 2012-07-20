subroutine da_ao_stats_sound (stats_unit, iv, re)

   !-------------------------------------------------------------
   ! Purpose: TBD   
   !-------------------------------------------------------------

   implicit none

   integer,        intent (in)    :: stats_unit    ! Output unit for stats.
   type (iv_type), intent (inout) :: iv            ! iv
   type (y_type),  intent (in)    :: re            ! A - O

   type (stats_sound_type) :: stats
   integer                 :: nu, nv, nt, nq
   integer                 :: n, k

   if (trace_use_dull) call da_trace_entry("da_ao_stats_sound")

   nu = 0
   nv = 0
   nt = 0
   nq = 0

   stats%maximum%u = maxmin_type (missing_r, 0, 0)
   stats%maximum%v = maxmin_type (missing_r, 0, 0)
   stats%maximum%t = maxmin_type (missing_r, 0, 0)
   stats%maximum%q = maxmin_type (missing_r, 0, 0)
   stats%minimum%u = maxmin_type(-missing_r, 0, 0)
   stats%minimum%v = maxmin_type(-missing_r, 0, 0)
   stats%minimum%t = maxmin_type(-missing_r, 0, 0)
   stats%minimum%q = maxmin_type(-missing_r, 0, 0)

   stats%average = residual_sound1_type(0.0, 0.0, 0.0, 0.0)
   stats%rms_err = stats%average

   do n=1, iv%info(sound)%nlocal
      do k=1, iv%info(sound)%levels(n)
         if (iv%info(sound)%proc_domain(1,n)) then
            call da_stats_calculate (n, k, iv%sound(n)%u(k)%qc, & 
               re%sound(n)%u(k), nu, &
               stats%minimum%u, stats%maximum%u, &
               stats%average%u, stats%rms_err%u)
            call da_stats_calculate (n, k, iv%sound(n)%v(k)%qc, & 
               re%sound(n)%v(k), nv, &
               stats%minimum%v, stats%maximum%v, &
               stats%average%v, stats%rms_err%v)
            call da_stats_calculate (n, k, iv%sound(n)%t(k)%qc, & 
               re%sound(n)%t(k), nt, &
               stats%minimum%t, stats%maximum%t, &
               stats%average%t, stats%rms_err%t)
            call da_stats_calculate (n, k, iv%sound(n)%q(k)%qc, & 
               re%sound(n)%q(k), nq, &
               stats%minimum%q, stats%maximum%q, &
               stats%average%q, stats%rms_err%q)
         end if
      end do
   end do

   ! Do inter-processor communication to gather statistics.
   call da_proc_sum_int (nu)
   call da_proc_sum_int (nv)
   call da_proc_sum_int (nt)
   call da_proc_sum_int (nq)
   iv%nstats(sound) = nu + nv + nt + nq

   call da_proc_stats_combine(stats%average%u, stats%rms_err%u, &
      stats%minimum%u%value, stats%maximum%u%value, &
      stats%minimum%u%n, stats%maximum%u%n, &
      stats%minimum%u%l, stats%maximum%u%l)
   call da_proc_stats_combine(stats%average%v, stats%rms_err%v, &
      stats%minimum%v%value, stats%maximum%v%value, &
      stats%minimum%v%n, stats%maximum%v%n, &
      stats%minimum%v%l, stats%maximum%v%l)
   call da_proc_stats_combine(stats%average%t, stats%rms_err%t, &
      stats%minimum%t%value, stats%maximum%t%value, &
      stats%minimum%t%n, stats%maximum%t%n, &
      stats%minimum%t%l, stats%maximum%t%l)
   call da_proc_stats_combine(stats%average%q, stats%rms_err%q, &
      stats%minimum%q%value, stats%maximum%q%value, &
      stats%minimum%q%n, stats%maximum%q%n, &
      stats%minimum%q%l, stats%maximum%q%l)

   if (rootproc) then
      if (nu /= 0 .or. nv /= 0 .or. nt /= 0 .or. nq /= 0) then
         write(unit=stats_unit, fmt='(/a/)') ' Diagnostics of AO for sound'
         call da_print_stats_sound(stats_unit, nu, nv, nt, nq, stats)
      end if
   end if

   if (trace_use_dull) call da_trace_exit("da_ao_stats_sound")

end subroutine da_ao_stats_sound

