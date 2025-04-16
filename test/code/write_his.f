      subroutine write_his(iloc)
c
c ----- write history file with global parameters
c       write tsfcc specially to preserve digits!
c
      use pars
      use fields
      use con_data
      use con_stats
c
      divgmax = 0.0
      do iz=1,nnz
         divgmax = amax1(divgmax, divz(iz))
      enddo
c
      ziavg = zi
      holtop = hol
      wt_min = wtsb(iloc,1)
      wt_le  = wtle(iloc,1)
      krec = krec + 1
      mid = nnz/4
      write(nhis1,6000) time,dt,utau,ziavg,amonin,holtop,
     +         (tsfcc(1)-t_ref),uwsfc,vwsfc,divgmax, wt_min, wt_le,
     +         ucfl, vcfl, wcfl, wtsfc(1),
     +         ups(mid),vps(mid),wps(mid),tps(mid,1),
     +         uwle(mid),uwsb(mid),uw_tot(mid),
     +         vwle(mid),vwsb(mid),vw_tot(mid),
     +         wtle(mid,1),wtsb(mid,1),wt_tot(mid,1),
     +         englez(mid),eavg(mid), wabs,
     +         u_10, v_10, cd_10, x_hurr, y_hurr,
     +         t_move, tau_x, tau_y, stokess, stokesw,
     +         cpou10, turb_la

 6000 format(5e17.8)
c
c -------------- write profile information
c
      call write_prof(nhisp,krec,isize,c_s%wwsb)
c
      return
      end
