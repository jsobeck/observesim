"""Master schedule module.

These routines allow the interpretation of the master schedule and the
calculation of basic astronomical parameters for each night within it.

The master schedule itself is kept as a Yanny file at:

 $OBSERVESIM_DIR/data/master_schedule.par

It contains a list of events labeled by local time and date, each of 
which either turns the survey on or turns it off. 

An example of how to use the master schedule is as follows:

import observesim.master as master

schedule = master.Master()
template = "mjd={mjd}, illumination={illumination}"
for mjd in schedule.mjds:
    night = master.Night(mjd=mjd)
    illumination = night.moon.illumination(mjd=mjd)
    print(template.format(mjd=mjd, illumination=illumination))

Dependencies:

 numpy
 scipy
 astropy
 pydl
 PyAstronomy

The reason for PyAstronomy is that it has Meeus-based routines for
astronomical calculations that are substantially faster than the more
accurate routines found in astropy.
"""

import os
import numpy as np
import scipy.optimize as optimize
import astropy.units as units
import astropy.time as atime
import pydl.pydlutils.yanny as yanny
import PyAstronomy.pyasl as pyasl
from observesim.moonphase import moonphase2
from observesim.sunpos2 import sunpos2


class MasterBase(object):
    """Master base class with generic internal utilities."""
    def __init__(self):
        return

    def _arrayify(self, quantity=None):
        """Cast quantity as ndarray of numpy.float64"""
        try:
            length = len(quantity)
        except TypeError:
            length = 1
        return np.zeros(length, dtype=np.float64) + quantity

    def _mjd2jd(self, mjd=None):
        """Convert MJD to JD"""
        return (self._arrayify(mjd) + np.float64(2400000.5))


class Observer(MasterBase):
    """Observer class to define different observatories.

    Parameters:
    ----------

    observatory : str
        Name of observatory to use (must be in observatory file)
        (default 'apo')

    observatoryfile : str
        Name of Yanny-format observatory file to read
        (default $OBSERVESIM_DIR/data/observatories.par)

    Attributes:
    ----------

    observatory : str
        Name of observatory

    latitude : numpy.float64
        Latitude of observatory

    longitude : numpy.float64
        Longitude (E of Greenwich) of observatory

    Methods:
    -------

    lst(mjd) : return LST in degrees for observer at given MJD (days)
"""
    def __init__(self, observatory='apo', observatoryfile=None):
        """Create Observer object"""
        self.observatory = observatory
        if(observatoryfile is None):
            observatoryfile = os.path.join(os.getenv('OBSERVESIM_DIR'),
                                           'data', 'observatories.par')
        self._file = observatoryfile
        self._data = yanny.yanny(self._file)
        observatories = np.array([observatory.decode()
                                  for observatory in
                                  self._data['OBSERVATORY']['observatory']])
        indx = np.where(observatories == self.observatory)[0][0]
        self.latitude = self._data['OBSERVATORY']['latitude'][indx]
        self.longitude = self._data['OBSERVATORY']['longitude'][indx]

    def lst(self, mjd=None):
        """Return LST (degrees) given MJD for observer"""
        mjds = self._arrayify(mjd)
        lst = (np.float64(15.) *
               pyasl.ct2lst(self._mjd2jd(mjds),
                            np.zeros(len(mjds)) + self.longitude))
        return (lst)


class Orb(MasterBase):
    """Orb class to handle spherical astronomy of objects in sky

    Parameters:
    ----------
    observer : Observer object
        Instance of Observer object (default an 'apo' Observer object)

    Methods:
    -------

    radec(mjd) : return (ra, dec) in deg J2000for Orb at MJD (days)
    hadec(mjd) : return (ha, dec) in deg of Orb for observer at MJD (days)
    altaz(mjd) : return (alt, az) in deg of Orb for observer at MJD (days)
"""
    def __init__(self, observer=None):
        """Create Orb object"""
        self.observer = observer
        if(self.observer is None):
            self.observer = Observer('apo')

    def radec(self, mjd=None):
        """Dummy function to report RA, Dec for generic Orb object"""
        return (np.float64(180.), np.float64(0.))

    def hadec(self, mjd=None):
        """Return (HA, dec) (degrees) of Orb given MJD for observer"""
        (ra, dec) = self.radec(self._arrayify(mjd))
        lst = self.observer.lst(mjd)
        ha = ((lst - ra + 360. + 180.) % 360.) - 180.
        return (ha, dec)

    def altaz(self, mjd=None):
        """Return (alt, az) (degrees) of Orb given MJD for observer"""
        # Need to make sure precession is dealh with correctly
        amjd = self._arrayify(mjd)
        (ha, dec) = self.hadec(amjd)
        (alt, az) = pyasl.hadec2altaz(ha, dec,
                                      np.zeros(len(amjd)) +
                                      self.observer.latitude)
        return (alt, az)


class Sun(Orb):
    """Sun class to track position of the Sun

    Parameters:
    ----------
    observer : Observer object
        Instance of Observer object (default an 'apo' Observer object)

    Methods:
    -------

    radec(mjd) : return (ra, dec) in deg J2000 for Sun at MJD (days)
    hadec(mjd) : return (ha, dec) in deg of Sun for observer at MJD (days)
    altaz(mjd) : return (alt, az) in deg of Sun for observer at MJD (days)
"""
    def radec(self, mjd=None):
        """Return (ra, dec) in deg J2000 for Sun at MJD (days)"""
        jd = self._mjd2jd(mjd=self._arrayify(mjd))
        (tmp_jd, ra, dec) = sunpos2(jd)
        return (ra, dec)


class Moon(Orb):
    """Moon class to track position and illumination of the Moon

    Parameters:
    ----------
    observer : Observer object
        Instance of Observer object (default an 'apo' Observer object)

    Methods:
    -------

    radec(mjd) : return (ra, dec) in deg J2000 for Moon at MJD (days)
    hadec(mjd) : return (ha, dec) in deg of Moon for observer at MJD (days)
    altaz(mjd) : return (alt, az) in deg of Moon for observer at MJD (days)
    illumination(mjd) : return illumination of Moon at MJD (days)
"""
    def radec(self, mjd=None):
        """Return (ra, dec) in deg J2000 for Moon at MJD (days)"""
        jd = self._mjd2jd(mjd=self._arrayify(mjd))
        ra, dec, dist, geolon, geolat = pyasl.moonpos(jd)
        return (ra, dec)

    def illumination(self, mjd=None):
        """Return Moon illumination at MJD (days)"""
        jd = self._mjd2jd(mjd=self._arrayify(mjd))
        return (moonphase2(jd))


class Night(MasterBase):
    """Night class to calculate observing parameters for night

    Parameters:
    ----------
    observer : Observer object
        Instance of Observer object (default an 'apo' Observer object)

    mjd : np.int32, int
        MJD of night (corresponds to MJD of morning)

    twilight : np.float64, float
        Twilight definition, in degrees above horizon (default -8)

    Attributes:
    ----------
    mjd_morning_twilight : np.float64
        MJD (days) of morning twilight

    mjd_evening_twilight : np.float64
        MJD (days) of evening twilight
"""

    def __init__(self, observer=None, mjd=None, twilight=-8.):
        self.mjd = np.int32(mjd)
        self.observer = observer
        self.twilight = twilight
        if(self.observer is None):
            self.observer = Observer('apo')
        self.sun = Sun(observer=self.observer)
        self.moon = Moon(observer=self.observer)
        self._calculate_evening_twilight()
        self._calculate_morning_twilight()

    def _twilight_function(self, mjd=None, twilight=-8.):
        (alt, az) = self.sun.altaz(mjd=mjd)
        return (alt - twilight)

    def _calculate_evening_twilight(self):
        noon_ish = (np.float64(self.mjd) -
                    self.observer.longitude / 15. / 24. - 0.5)
        midnight_ish = noon_ish + 0.5
        twi = optimize.brenth(self._twilight_function,
                              noon_ish, midnight_ish,
                              args=(self.twilight))
        self.mjd_evening_twilight = np.float64(twi)

    def _calculate_morning_twilight(self):
        midnight_ish = (np.float64(self.mjd) -
                        self.observer.longitude / 15. / 24.)
        nextnoon_ish = midnight_ish + 0.5
        twi = optimize.brenth(self._twilight_function,
                              midnight_ish, nextnoon_ish,
                              args=(self.twilight))
        self.mjd_morning_twilight = np.float64(twi)


class Master(MasterBase):
    """Master class to interpret master schedule

    Parameters:
    ----------
    schedulefile : str
        schedule file to use; default $OBSERVESIM_DIR/data/master_schedule.par

    Attributes:
    ----------
    start : np.int32
        MJD (days) of first night of survey
    end : np.int32
        MJD (days) of last night of survey
    mjds : ndarray of np.int32
        MJDs (days) when survey is potentially active
    events : ndarray of numpy.str_
        names of events of note
    event_dates : ndarray of numpy.str_
        list of dates in ISO format for events of note
    event_times : ndarray of numpy.str_
        list of times in ISO format for events of note
    event_mjd : ndarray of numpy.float64
        MJDs (days) of events of note

    Methods:
    -------
    on() : is the survey on
"""
    def __init__(self, schedulefile=None):
        """Create Master object for schedule"""
        if(schedulefile is None):
            schedulefile = os.path.join(os.getenv('OBSERVESIM_DIR'),
                                        'data', 'master_schedule.par')
        self._schedulefile = schedulefile
        self.schedule = yanny.yanny(self._schedulefile)
        self._validate()
        self.event_dates = np.array([date.decode() for date
                                     in self.schedule['SCHEDULE']['date']])
        self.event_times = np.array([time.decode() for time
                                     in self.schedule['SCHEDULE']['time']])
        self.event_mjds = self._dateandtime2mjd()
        self.events = np.array([event.decode() for event
                                in self.schedule['SCHEDULE']['event']])
        self.start = self._start()
        self.end = self._end()
        self.mjds = self._mjds()
        return

    def _dateandtime2mjd(self):
        isotimes = ["{date} {time}".format(date=date, time=time) for date, time
                    in zip(self.event_dates, self.event_times)]
        times = atime.Time(isotimes, format='iso', scale='tai')
        times = times + np.int32(self.schedule['to_tai']) * units.hour
        return(times.mjd)

    def _validate(self):
        # should make sure:
        #  one start (first event)
        #  one end (last event)
        #  start MJD is a daytime time 
        #  START_SURVEY is "on" 
        #  END_SURVEY is "off" 
        return

    def on(self, mjd=None):
        if(mjd < self.event_mjds[0]):
            return('off', self.event_mjds[0])
        if(mjd >= self.event_mjds[-1]):
            return('off', mjd + 1.)
        # Assumes there is only one
        indx = np.where((mjd >= self.event_mjds[0:-1]) &
                        (mjd < self.event_mjds[1:]))[0][0]
        return(self.schedule[self.events[indx]],
               self.event_mjds[indx + 1])

    def end_mjd(self):
        return(self.event_mjds[-1])

    def _start(self):
        # Assumes there is only one
        indx = np.where(self.events == 'START_SURVEY')[0][0]
        # Assumes START_SURVEY turns state on
        return(np.int32(np.floor(self.event_mjds[indx])))

    def _end(self):
        # Assumes there is only one
        indx = np.where(self.events == 'END_SURVEY')[0][0]
        # Assumes END_SURVEY turns state off
        return(np.int32(np.ceil(self.event_mjds[indx])))

    def _mjds(self):
        nmjd = self.end - self.start + 1
        mjds = self.start + np.arange(nmjd, dtype=np.int32)
        keep = np.zeros(nmjd, dtype=np.int32)
        for indx in np.arange(len(self.events) - 1):
            this_event = self.events[indx]
            if(self.schedule[this_event] == 'on'):
                keep_start = np.int32(np.floor(self.event_mjds[indx]))
                keep_end = np.int32(np.ceil(self.event_mjds[indx + 1]))
                ikeep = np.where((mjds >= keep_start) &
                                 (mjds <= keep_end))[0]
                keep[ikeep] = 1
        ikeep = np.where(keep)[0]
        return(mjds[ikeep])
