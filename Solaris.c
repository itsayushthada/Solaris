/* For each line if it is asking for more than one input  then give the values seprated by space*/
#include <stdio.h>
#include <math.h>
#define UTC 8
/*
UTC Values For Various Countries:
India: 5.5
Egypt: 2
Moscow: 3
Berlin: 2
Rio De Janeiro: -3
Columbia: -5
Chicago: -6
Quebec: -4
Vancouver: -7
China: 8
Japan: 9
Norway: 2
Madagascar: 3
Iraq: 3
Denmark: -3
Lisbon: 1
Ireland: 1
Marshall Islands: 12
Brisbane: 10
Anchorage: -9
*/

/* A macro to compute the number of days elapsed since 0.0 Jan 2000 */
/* (which is equivalent to 1999 Dec 31, 0h UT)                      */

#define days_since_2000_Jan_0(y,m,d) \
    (367L*(y)-((7*((y)+(((m)+9)/12)))/4)+((275*(m))/9)+(d)-730530L)

/* Conversion factors radians and degrees */

#ifndef PI
#define PI        3.1415926535897932384
#endif

#define RadianToDegree     ( 180.0 / PI )
#define DegreeToRadian    ( PI / 180.0 )

/* The trigonometric functions in degrees */

#define sind(x)  sin((x)*DegreeToRadian)
#define cosd(x)  cos((x)*DegreeToRadian)
#define tand(x)  tan((x)*DegreeToRadian)

#define atand(x)    (RadianToDegree*atan(x))
#define asindh(x)    (RadianToDegree*asin(x))
#define acosd(x)    (RadianToDegree*acos(x))
#define atan2d(y,x) (RadianToDegree*atan2(y,x))

/* Following are some macros around the "workhorse" function __daylen__ */
/* They mainly fill in the desired values for the reference altitude    */
/* below the horizon, and also selects whether this altitude should     */
/* refer to the Sun's center or its upper limb.                         */

/* This macro computes the length of the day, from sunrise to sunset.   */
/* Sunrise/Sunset is considered to occur when the Sun's upper limb is   */
/* 35 arc minutes below the horizon (this accounts for the refraction   */
/* of the Earth's atmosphere).                                          */
/* Here the horizon is taken as 0 degrees.                               */
#define day_length(year,month,day,lon,lat)  \
        __daylen__( year, month, day, lon, lat, -35.0/60.0, 1 )

/* This macro computes the length of the day, including civil twilight. */
/* Civil twilight starts/ends when the Sun's center is 6 degrees below  */
/* the horizon.                                                         */
/* Here the horizon is taken as 0 degrees.                               */
#define day_civil_twilight_length(year,month,day,lon,lat)  \
        __daylen__( year, month, day, lon, lat, -6.0, 0 )

/* This macro computes the length of the day, incl. nautical twilight.  */
/* Nautical twilight starts/ends when the Sun's center is 12 degrees    */
/* below the horizon.                                                   */
/* Here the horizon is taken as 0 degrees.                               */
#define day_nautical_twilight_length(year,month,day,lon,lat)  \
        __daylen__( year, month, day, lon, lat, -12.0, 0 )

/* This macro computes the length of the day, incl. astronomical twilight. */
/* Astronomical twilight starts/ends when the Sun's center is 18 degrees   */
/* below the horizon.                                                      */
/* Here the horizon is taken as 0 degrees.                               */
#define day_astronomical_twilight_length(year,month,day,lon,lat)  \
        __daylen__( year, month, day, lon, lat, -18.0, 0 )

/* This macro computes times for sunrise/sunset.                      */
/* Sunrise/set is considered to occur when the Sun's upper limb is    */
/* 35 arc minutes below the horizon (this accounts for the refraction */
/* of the Earth's atmosphere).                                        */
/* Here the horizon is taken as 0 degrees.                               */
#define sun_rise_set(year,month,day,lon,lat,rise,set)  \
        __sunriset__( year, month, day, lon, lat, -35.0/60.0, 1, rise, set )

/* This macro computes the start and end times of civil twilight.       */
/* Civil twilight starts/ends when the Sun's center is 6 degrees below  */
/* the horizon.                                                         */
/* Here the horizon is taken as 0 degrees.                               */
#define civil_twilight(year,month,day,lon,lat,start,end)  \
        __sunriset__( year, month, day, lon, lat, -6.0, 0, start, end )

/* This macro computes the start and end times of nautical twilight.    */
/* Nautical twilight starts/ends when the Sun's center is 12 degrees    */
/* below the horizon.                                                   */
/* Here the horizon is taken as 0 degrees.                               */
#define nautical_twilight(year,month,day,lon,lat,start,end)  \
        __sunriset__( year, month, day, lon, lat, -12.0, 0, start, end )

/* This macro computes the start and end times of astronomical twilight.   */
/* Astronomical twilight starts/ends when the Sun's center is 18 degrees   */
/* below the horizon.                                                      */
/* Here the horizon is taken as 0 degrees.                               */
#define astronomical_twilight(year,month,day,lon,lat,start,end)  \
        __sunriset__( year, month, day, lon, lat, -18.0, 0, start, end )

/* Declaring function signatures */

double __daylen__( int year, int month, int day, double lon, double lat,double altit, int upper_limb );
int __sunriset__( int year, int month, int day, double lon, double lat,double altit, int upper_limb, double *rise, double *set );
void sunpos( double d, double *lon, double *r );
void sun_RA_dec( double d, double *RA, double *dec, double *r );
double revolution( double x );
double rev180( double x );
double GMST0( double d );
void convert_am_pm(float std_time,double value);

int main()
{
    int year,month,day;
    double lon, lat;
    double daylen, civlen, nautlen, astrlen;
    double rise, set, civ_start, civ_end, naut_start, naut_end, astr_start, astr_end;
    int    rs, civ, naut, astr;
    char 	 buf[80];

    printf( "Longitude(+E) and Latitude(+N) :\t" );
    fgets(buf, 80, stdin);
    sscanf(buf, "%lf %lf", &lon, &lat );

    while(1)
    {
        printf( "Enter the Required Date (DD-MM-YYYY):\t" );
        fgets(buf, 80, stdin);
        sscanf(buf, "%d %d %d", &day , &month, &year);

        daylen  = day_length(year,month,day,lon,lat);
        civlen  = day_civil_twilight_length(year,month,day,lon,lat);
        nautlen = day_nautical_twilight_length(year,month,day,lon,lat);
        astrlen = day_astronomical_twilight_length(year,month,day,lon,lat);

        printf( "\n\nDay length:\t\t\t\t");
        convert_am_pm(0,daylen);
        printf( "\nWith Civil Twilight:\t\t\t");
        convert_am_pm(0,civlen );
        printf( "\nWith Nautical Twilight:\t\t\t");
        convert_am_pm(0,nautlen)  ;
        printf( "\nWith Astronomical Twilight:\t\t");
        convert_am_pm(0,astrlen);
        printf( "\n\nLength of Civil Twilight:\t\t");
        convert_am_pm(0,(civlen-daylen)/2.0);
        printf( "\nLength of Nautical Twilight:\t\t");
        convert_am_pm(0,(nautlen-daylen)/2.0);
        printf( "\nLength of Astronomical Twilight:\t");
        convert_am_pm(0,(astrlen-daylen)/2.0);

        rs   = sun_rise_set( year, month, day, lon, lat,&rise, &set );
        civ  = civil_twilight( year, month, day, lon, lat,&civ_start, &civ_end );
        naut = nautical_twilight( year, month, day, lon, lat,&naut_start, &naut_end );
        astr = astronomical_twilight( year, month, day, lon, lat, &astr_start, &astr_end );

        printf( "\n\nSun at South: ");
        convert_am_pm(UTC,(rise+set)/2.0) ;
        switch( rs )
        {
            case 0:
                    printf( "\n\nSun Rises:\t\t\t\t");
					convert_am_pm(UTC,rise);
					printf( "\nSun Sets:\t\t\t\t");
					convert_am_pm(UTC,set);
                    break;
            case +1:
                    printf( "\nSun above Horizon: " );
                    break;
            case -1:
                    printf( "\nSun below Horizon: " );
                    break;
        }
        switch( civ )
        {
            case 0:
                    printf( "\n\nCivil Twilight Starts:\t\t\t");
                    convert_am_pm(UTC,civ_start);
                    printf("\nCivil Twilight Ends:\t\t\t");
					convert_am_pm(UTC,civ_end);
                    break;
            case +1:
                    printf( "\nNever darker than Civil twilight: " );
                    break;
            case -1:
                    printf( "\nNever as bright as Civil twilight: " );
                    break;
        }
        switch( naut )
        {
            case 0:
                    printf( "\n\nNautical Twilight Starts:\t\t");
                    convert_am_pm(UTC,naut_start);
                    printf( "\nNautical Twilight Ends:\t\t\t");
					convert_am_pm(UTC,naut_end);
                    break;
            case +1:
                    printf( "\nNever darker than nautical twilight: " );
                    break;
            case -1:
                    printf( "\nNever as bright as nautical twilight: " );
                    break;
        }
        switch( astr )
        {
            case 0:
                    printf( "\n\nAstronomical Twilight Starts:\t\t");
                    convert_am_pm(UTC,astr_start);
                    printf( "\nAstronomical Twilight Ends\t\t");
					convert_am_pm(UTC,astr_end);
                    break;
            case +1:
                    printf( "\nNever darker than Astronomical Twilight: " );
                    break;
            case -1:
                    printf( "\nNever as bright as Astronomical Twilight: " );
                    break;
        }
      printf("\n\n");
      return 0;
      }
}

/* The "workhorse" function for Sunrise/Sunset timings */

int __sunriset__( int year, int month, int day, double lon, double lat, double altit, int upper_limb, double *trise, double *tset )
/***************************************************************************/
/* Note: year,month,date = calendar date, from 1801 to 2099 only.     */
/*       Take Eastern longitude positive & Western longitude negative */
/*       Take Northern latitude positive & Southern latitude negative */
/*       The longitude value IS critical in this function!            */
/*       altit = the altitude which the Sun should cross              */
/*               Set to -35/60 degrees for sunrise/sunset, -6 degrees */
/*               for civil, -12 degrees for nautical and -18          */
/*               degrees for astronomical twilight respectively.      */
/*         upper_limb: non-zero -> upper limb, zero -> center         */
/*               Set to non-zero (e.g. 1) when computing sunrise/sunset */
/*               times, and to zero when computing start/end of       */
/*               twilight.                                            */
/*        *rise = where to store the rise time                        */
/*        *set  = where to store the set  time                        */
/*                Both times are relative to the specified altitude,  */
/*                and thus this function can be used to compute       */
/*                various twilight times, as well as rise/set times   */
/* Return value:  0 = sun rises/sets this day, times stored at        */
/*                    *trise and *tset.                               */
/*               +1 = sun above the specified "horizon" 24 hours.     */
/*                    *trise set to time when the sun is in the south,    */
/*                    minus 12 hours while *tset is set to the south  */
/*                    time plus 12 hours. "Day" length = 24 hours     */
/*               -1 = sun is below the specified "horizon" 24 hours   */
/*                    "Day" length = 0 hours, *trise and *tset are    */
/*                    both set to the time when the sun is in the south.  */
/*                                                                    */
/**********************************************************************/
{
      double  d,  /* Days since 2000 Jan 0.0 (value is negative before 2000 Jan 0.0) */
      sr,         /* Solar Distance, Astronomical Units */
      sRA,        /* Sun's Right Ascension */
      sdec,       /* Sun's Declination */
      sradius,    /* Sun's Apparent Radius */
      t,          /* Diurnal Arc */
      tsouth,     /* Time When Sun Is In The South */
      sidtime;    /* Local Sidereal Time */

      int rc = 0; /* Return cde from function - usually 0 */

      /* Compute d of 12h local mean solar time */
      d = days_since_2000_Jan_0(year,month,day) + 0.5 - lon/360.0;

      /* Compute the local sidereal time of this moment */
      sidtime = revolution( GMST0(d) + 180.0 + lon );

      /* Compute Sun's RA, Decl and distance at this moment */
      sun_RA_dec( d, &sRA, &sdec, &sr );

      /* Compute time when Sun is at south - in hours UT */
      tsouth = 12.0 - rev180(sidtime - sRA)/15.0;

      /* Compute the Sun's apparent radius in degrees */
      sradius = 0.2666 / sr;

      /* Do correction to upper limb, if necessary */
      if ( upper_limb )
            altit -= sradius;

      /* Compute the diurnal arc that the Sun traverses to reach */
      /* the specified altitude altitude: */
      {
            double cost;
            cost = ( sind(altit) - sind(lat) * sind(sdec) ) /
                  ( cosd(lat) * cosd(sdec) );
            if ( cost >= 1.0 )
                  rc = -1, t = 0.0;       /* Sun always below altitude */
            else if ( cost <= -1.0 )
                  rc = +1, t = 12.0;      /* Sun always above altitude */
            else
                  t = acosd(cost)/15.0;   /* The diurnal arc, hours */
      }
      /* Store rise and set times - in hours UT */
      *trise = tsouth - t;
      *tset  = tsouth + t;

      return rc;
}
/* __sunriset__ */

/* The "workhorse" function */
double __daylen__( int year, int month, int day, double lon, double lat,double altit, int upper_limb )
/**********************************************************************/
/* Note: year,month,date = calendar date, from 1801 to 2099 only.     */
/*       Take Eastern longitude positive & Western longitude negative */
/*       Take Northern latitude positive & Southern latitude negative */
/*       The longitude value is not critical. Set it to the correct   */
/*       longitude if you're picky, otherwise set to to, say, 0.0     */
/*       The latitude however IS critical - be sure to get it correct */
/*       altit = the altitude which the Sun should cross              */
/*               Set to -35/60 degrees for rise/set, -6 degrees       */
/*               for civil, -12 degrees for nautical and -18          */
/*               degrees for astronomical twilight respecrively.      */
/*         upper_limb: non-zero -> upper limb, zero -> center         */
/*               Set to non-zero (e.g. 1) when computing day length   */
/*               and to zero when computing day & twilight length combined.      */
/**********************************************************************/
{
      double  d,  /* Days since 2000 Jan 0.0 (value is negative before 2000 Jan 0.0) */
      obl_ecl,    /* Obliquity (Inclination) Of Earth's Axis */
      sr,         /* Solar Distance, Astronomical Units */
      slon,       /* True Solar Longitude */
      sin_sdecl,  /* Sine Of Sun's Declination */
      cos_sdecl,  /* Cosine Of Sun's Declination */
      sradius,    /* Sun's Apparent Radius */
      t;          /* Diurnal Arc */

      /* Compute d of 12h local mean solar time */
      d = days_since_2000_Jan_0(year,month,day) + 0.5 - lon/360.0;

      /* Compute obliquity of ecliptic (inclination of Earth's axis) */
      obl_ecl = 23.4393 - 3.563E-7 * d;

      /* Compute Sun's ecliptic longitude and distance */
      sunpos( d, &slon, &sr );

      /* Compute sine and cosine of Sun's declination */
      sin_sdecl = sind(obl_ecl) * sind(slon);
      cos_sdecl = sqrt( 1.0 - sin_sdecl * sin_sdecl );

      /* Compute the Sun's apparent radius, degrees */
      sradius = 0.2666 / sr;

      /* Do correction to upper limb, if necessary */
      if ( upper_limb )
            altit -= sradius;

      /* Compute the diurnal arc that the Sun traverses to reach */
      /* the specified altitude altit: */
      {
            double cost;
            cost = ( sind(altit) - sind(lat) * sin_sdecl ) /
                  ( cosd(lat) * cos_sdecl );
            if ( cost >= 1.0 )
                  t = 0.0;                      /* Sun always below altit */
            else if ( cost <= -1.0 )
                  t = 24.0;                     /* Sun always above altit */
            else  t = (2.0/15.0) * acosd(cost); /* The diurnal arc, hours */
      }
      return t;
}
/* __daylen__ */

/* This function computes the Sun's position at any given instant */

void sunpos( double d, double *lon, double *r )
/******************************************************/
/* Computes the Sun's ecliptic longitude and distance */
/* at an instant given in d, number of days since     */
/* 2000 Jan 0.0.  The Sun's ecliptic latitude is not  */
/* computed, since it's always very near 0.           */
/******************************************************/
{
      double M,         /* Mean Anomaly Of The Sun */
             w,         /* Mean Longitude Of Perihelion */
                        /* Note: Sun's Mean Longitude = M + w */
             e,         /* Eccentricity Of Earth's Orbit */
             E,         /* Eccentric Anomaly */
             x, y,      /* X, Y Coordinates In Orbit */
             v;         /* True Anomaly */

      /* Compute mean elements */
      M = revolution( 356.0470 + 0.9856002585 * d );
      w = 282.9404 + 4.70935E-5 * d;
      e = 0.016709 - 1.151E-9 * d;

      /* Compute true longitude and radius vector */
      E = M + e * RadianToDegree * sind(M) * ( 1.0 + e * cosd(M) );
            x = cosd(E) - e;
      y = sqrt( 1.0 - e*e ) * sind(E);
      *r = sqrt( x*x + y*y );              /* Solar distance */
      v = atan2d( y, x );                  /* True anomaly */
      *lon = v + w;                        /* True solar longitude */
      if ( *lon >= 360.0 )
            *lon -= 360.0;                   /* Make it 0..360 degrees */
}

void sun_RA_dec( double d, double *RA, double *dec, double *r )
/******************************************************/
/* Computes the Sun's equatorial coordinates RA, Decl */
/* and also its distance, at an instant given in d,   */
/* the number of days since 2000 Jan 0.0.             */
/******************************************************/
{
      double lon, obl_ecl, x, y, z;

      /* Compute Sun's ecliptical coordinates */
      sunpos( d, &lon, r );

      /* Compute ecliptic rectangular coordinates (z=0) */
      x = *r * cosd(lon);
      y = *r * sind(lon);

      /* Compute obliquity of ecliptic (inclination of Earth's axis) */
      obl_ecl = 23.4393 - 3.563E-7 * d;

      /* Convert to equatorial rectangular coordinates - x is unchanged */
      z = y * sind(obl_ecl);
      y = y * cosd(obl_ecl);

      /* Convert to spherical coordinates */
      *RA = atan2d( y, x );
      *dec = atan2d( z, sqrt(x*x + y*y) );
}
/* sun_RA_dec */

/******************************************************************/
/* This function reduces any angle to within the first revolution */
/* by subtracting or adding even multiples of 360.0 until the     */
/* result is >= 0.0 and < 360.0                                   */
/******************************************************************/

#define INV360    ( 1.0 / 360.0 )

double revolution( double x )
/*****************************************/
/* Reduce angle to within 0..360 degrees */
/*****************************************/
{
      return( x - 360.0 * floor( x * INV360 ) );
}  /* revolution */

double rev180( double x )
/*********************************************/
/* Reduce angle to within +180..+180 degrees */
/*********************************************/
{
      return( x - 360.0 * floor( x * INV360 + 0.5 ) );
}  /* revolution */

/*******************************************************************/
/* This function computes GMST0, the Greenwich Mean Sidereal Time  */
/* at 0h UT (i.e. the sidereal time at the Greenwhich meridian at  */
/* 0h UT).  GMST is then the sidereal time at Greenwich at any     */
/* time of the day.  I've generalized GMST0 as well, and define it */
/* as:  GMST0 = GMST - UT  --  this allows GMST0 to be computed at */
/* other times than 0h UT as well.  While this sounds somewhat     */
/* contradictory, it is very practical:  instead of computing      */
/* GMST like:                                                      */
/*                                                                 */
/*  GMST = (GMST0) + UT * (366.2422/365.2422)                      */
/*                                                                 */
/* where (GMST0) is the GMST last time UT was 0 hours, one simply  */
/* computes:                                                       */
/*                                                                 */
/*  GMST = GM
ST0 + UT                                              */
/*                                                                 */
/* where GMST0 is the GMST "at 0h UT" but at the current moment!   */
/* Defined in this way, GMST0 will increase with about 4 min a     */
/* day.  It also happens that GMST0 (in degrees, 1 hr = 15 degr)   */
/* is equal to the Sun's mean longitude plus/minus 180 degrees!    */
/* (if we neglect aberration, which amounts to 20 seconds of arc   */
/* or 1.33 seconds of time)                                        */
/*                                                                 */
/*******************************************************************/

double GMST0( double d )
{
      double sidtim0;
      /* Sidtime at 0h UT = L (Sun's mean longitude) + 180.0 degr  */
      /* L = M + w, as defined in sunpos().                        */
      sidtim0 = revolution( ( 180.0 + 356.0470 + 282.9404 ) +
                          ( 0.9856002585 + 4.70935E-5 ) * d );
      return sidtim0;
}  /* GMST0 */

void convert_am_pm(float std_time,double value)
{
	int n,hr,min,sec;
	n=(std_time+value)*60*60;
	if(n>3600){
		min = n/60;
		sec = n%60;
		hr = min/60;
		min = min%60;
		printf("%d:%d:%d",hr,min,sec);
	}
	else{
		min = n/60;
		sec = n%60;
		printf("0:%d:%d",min,sec);
	}
}
