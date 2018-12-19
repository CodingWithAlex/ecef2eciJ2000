/*
    This class transforms a vector from the earth-centered earth-fixed (ecef/itrf) frame, 
    to the earth-centered intertial (eci) mean-equator mean-equinox (j2000) frame, using
    full precession, nutation, rotation, polar motion intermediaries.

    Ref: Vallado, Fundamentals of Astrodynamics and Applications
    
    This code is largely based on (Matlab) code provided by Vallado: https://celestrak.com/software/vallado-sw.php
*/
package ca.gc.asccsa.neossat.util;

import ca.gc.asccsa.neossat.FITSProcessor.FitsProcessorArgs;
import ca.gc.asccsa.neossat.dataStructures.TIME42;
import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.io.FileUtils;

/**
 *
 * @author SHasaj (Semi Hasaj, CSA student)
 */
public class ecef2ecij2000 {
    // **** EOP.dat **** //
    static String EOP_PATH;
    static double dUT1 = 0;                       // diff between UTC and UT1-- UT1 = UTC + dUT1, [sec]
    static double xp =  0;                                  // polar motion coefficient, x, [arcsec] 
    static double yp = 	0;                                  // polar motion coefficient, y, [arcsec]
    // ***************** //
    
    static String CONFIG_PATH;
    public static double lod = 0;                              // excess length of day, [sec]
    public static double ddpsi = 0;                               // delta psi correction to gcrf, [marcsec]
    public static double ddeps = 0;                               // delta eps correction to gcrf, [marcsec]

    static String LEAPSEC_PATH;
    public static int dAT = 0;                                // diff between UTC and TAI-- TAI = UTC + dAT, [sec]
    
    // **** nut80.dat **** //
    static String NUT80_PATH;
    static double[][] iar80 = new double[106][5];            // integers for fk5 1980
    static double[][] rar80 = new double[106][4];            // reals for fk5 1980, [rad]
    // ***************** //
    
    public static boolean isIAU80_loaded = false;
    public static boolean isEOP_loaded = false;
    public static boolean isFitsConfig_loaded = false;
    public static double MJD = 0;
    public int eqeterms = 0;                        // terms for sideral AST calc. [0 or 2]
    public final double PI = Math.PI; 
        
    public double[][] prec;                         // rotation for precession
    public double[][] nut;                          // rotation for nutation
    public double[][] st;                           // rotation for sidereal time
    public double[][] stdot = new double[3][3];     // transformation for pef-tod rate
    public double[][] pm;                           // rotation for polar motion
    public double deltapsi;                         // nutation angle
    public double meaneps;                          // mean obliquity of the ecliptic
    public double omega;                            // longitude of ascending node of moon
    
    // ECI OUTPUTS
    public double[] r_eci;
    public double[] r_pef;
    public double[] v_eci;
    public double[] v_pef;
    
    public ecef2ecij2000(double[] r_ecef, double[] v_ecef,  String timeUTC) throws ParseException {
        // time format: "yyyy-DDD-HH:mm:ss.SSS" UTC
        this(r_ecef, v_ecef, new TIME42(timeUTC, new SimpleDateFormat("yyyy-DDD-HH:mm:ss.SSS")));
    }
        
    public ecef2ecij2000(double[] r_ecef, double[] v_ecef,  TIME42 time42) {

        // time42 is in UTC
        double jdut1 = time42.getJD_AsDouble();         // julian date of UTC
        
        double mjd_current = jdut1-2400000.5;
        // if we've advanced a full julian day, we need to grab a new set of data from eop.dat
        if (Math.abs(MJD - mjd_current) >= 1) {
            // revert the eop.data flag so that initialize() will know to re-process the eop.dat file with a new MJD to get more up-to-date values
            MJD = mjd_current;
            isEOP_loaded = false;
        }
        
        // set constants that are read from files;
        initialize();
        
        jdut1 += (dUT1) / 86400.0;                      // julian date of UT1, accounting for sub-millisecond timing of dUT1
        dAT = getLeapSecond(LEAPSEC_PATH, jdut1);
        time42.addMillis((long) (dAT*1000 + 32184));    // Terrestrial Time (TT) = dTAI + dGPS + 32.184 = 19 + dGPS + 32.184 = dAT + 32.184
        double JD = time42.getJD_AsDouble();            // JD of TT
        double ttt = (JD - 2451545.0) / 36525;          // julian centuries of TT, from J2000 epoche
        
        prec = precession(ttt);
        nut = nutation(ttt, ddpsi, ddeps);
        st = sidereal(jdut1);
        pm = polarMotion(xp, yp, ttt);
        
        double thetasa = 7.29211514670698e-05 * (1.0  - lod/86400.0 );
        double[] omegaearth = {0, 0, thetasa}; 
        r_pef = matrixMult(pm, r_ecef);
        double[] r_eci_temp = matrixMult(st, r_pef);
        r_eci_temp = matrixMult(nut, r_eci_temp);
        r_eci_temp = matrixMult(prec, r_eci_temp);
        r_eci = r_eci_temp;
        
        v_pef = matrixMult(pm, v_ecef);
        double[] v_eci_temp = matrixMult( st, VectorAdd(v_pef, crossProduct(omegaearth, r_pef)) );
        v_eci_temp = matrixMult(nut, v_eci_temp);
        v_eci_temp = matrixMult(prec, v_eci_temp);
        v_eci = v_eci_temp;
    }
    
    public static void printMatrix(double[][] m) {
        for (int i = 0; i < m.length; i++) {
            double[] ds = m[i];
            for (int j = 0; j < ds.length; j++) {
                double d = ds[j];
                System.out.print(d + "  ");
            }
            System.out.println();
        }
    }
    
    public void initialize() {
        NUT80_PATH = FitsProcessorArgs.getInstance().getNut80Path();
        EOP_PATH = FitsProcessorArgs.getInstance().getEopPath();
        LEAPSEC_PATH = FitsProcessorArgs.getInstance().getLeapSecondPath();
        CONFIG_PATH = FitsProcessorArgs.getInstance().getFitsConfigPath();
        if (!isIAU80_loaded) {
            try {
                File f = new File(NUT80_PATH);
                if (!f.exists()) {
                    System.out.println("no nut80.dat file found at path: " + NUT80_PATH);
                    return;
                }
                List<String> ls = getFileContent(f);
                for (int i = 0; i < ls.size(); i++) {
                    String string = ls.get(i).trim();
                    String[] s = string.split("(\\s)+");
                    for (int j = 0; j < s.length-1; j++) {
                        if (j < 5) { // first 5 columns contain data for iar80, next 5 for rar80
                            iar80[i][j] = Double.valueOf(s[j]);
                        }
                        else {
                            rar80[i][j-5] = arcsec2rad( Double.valueOf(s[j]) / 10000 ); // 0.0001" to rad
                        }
                    }
                }
                isIAU80_loaded = true;
            } catch (IOException ex) {
                Logger.getLogger(ecef2ecij2000.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
        if(!isEOP_loaded) {
            try {
                File f2 = new File(EOP_PATH);
                if (!f2.exists()) {
                    System.out.println("no " +FitsProcessorArgs.DEFAULT_EOP_NAME+ " file found at path: " + EOP_PATH);
                    return;
                }
                List<String> lst = getFileContent(f2);
                double last_lod_MJD = 0;
                for (String row : lst) {
//                    String str = lst.get(i).substring(7).trim(); // first 6 characters represent yyMMdd, with no consisting spacing, so we just skip them altogether
//                    String[] s = str.split("(\\s)+");
                    double mjd = Double.valueOf( row.substring(6, 16) ); // first column is MJD
                    // lod generally is not available until a few days after a finals.daily file is produced, so we used the closest provided lod value
                    String lodString = row.substring(79, 87).trim();
                    if (!"".equals(lodString)) {
                        lod = Double.valueOf(lodString) / 1000.0; // msec -> sec
                        last_lod_MJD = mjd;
                    }
                    if (Math.abs(mjd - MJD) <= 0.5) {
                        // this represents the closest MJD to the MJD of the current coordinates
                        xp = Double.valueOf(row.substring(18, 28));
                        yp = Double.valueOf(row.substring(37, 47));
                        dUT1 = Double.valueOf(row.substring(58, 69));
                        ddpsi = Double.valueOf(row.substring(98, 107)) / 1000.0; // milliarcsec -> arcsec
                        ddeps = Double.valueOf(row.substring(118, 126)) / 1000.0; // milliarcsec -> arcsec
                        
                        // accept the last recorded lod if it less than a month (30 days) old. Else, set it equal 0
                        if ( Math.abs(last_lod_MJD - mjd) > 30 ) {
                            lod = 0;
                        }
                        break;
                    }
                }
                isEOP_loaded = true;
            } catch (Exception e) {
                System.err.println("Error extracting finals.daily file");
            }
        }
        
    }
    
    public static int getLeapSecond(String leap_second_file_path, double givenJD) {
        try {
            int leapSec = 0;
            File f = new File(leap_second_file_path);
            if (!f.exists()) {
                System.err.println("no LeapSecond.dat file found at path: " + leap_second_file_path);
                return 37;  // default latest value as of Dec 18, 2018
            }
            List<String> ls = getFileContent(f);
            for (String l : ls) {
                String[] columns = l.trim().split("(\\s)+");
                if (columns.length < 7) {
                    continue;
                }
                String julianDate = columns[3].trim();
                double jd = Double.valueOf(julianDate);
                if (givenJD < jd) {
                    // compare the given JD with that in the current LeapSecond.dat row.
                    // if we see a JD in this file that is a time after the given JD, 
                    // return the previously recorded leap second value, as this is the correct
                    // leap second for the given JD.
                    return leapSec;
                }
                String leapsec = columns[4].trim();
                leapSec = Double.valueOf(leapsec).intValue();
            }
            return leapSec; // reached when the given JD is after the last recorded Leap Second update.
        } catch (Exception e) {
            return 37;  // default latest value as of Dec 18, 2018
        }
    }
    
    public static void initiateConfig() {
        try {
                    File f2 = new File(CONFIG_PATH);
                    if (!f2.exists()) {
                        System.out.println("no fitsconfig.txt file found at path: " + CONFIG_PATH);
                        dAT = 37; // default value as of Nov. 20, 2018
                        return;
                    }
                    List<String> lst = getFileContent(f2);
                    for (int i = 0; i < lst.size(); i++) {
                        String str = lst.get(i).trim();
                        String[] s = str.split("(\\s)+");
                        String name = s[0].trim();
                        if (name.equals("dAT")) {
                            dAT = Integer.parseInt(s[2].trim() );
                        }
                    }
                    isFitsConfig_loaded = true;
                } catch (Exception e) {
                    System.err.println("Error extracting fitsconfig.txt file");
                    e.printStackTrace();
                    dAT = 37; // default value as of Nov. 20, 2018
                }
    }
    
    public double [][] polarMotion(double xp, double yp, double ttt) {
        double[][] pm = new double[3][3];
        xp = arcsec2rad(xp);
        yp = arcsec2rad(yp);
        double cosxp = Math.cos(xp);
        double sinxp = Math.sin(xp);
        double cosyp = Math.cos(yp);
        double sinyp = Math.sin(yp);
        
            pm[0][0] =  cosxp;
            pm[0][1] =  0.0;
            pm[0][2] = -sinxp;
            pm[1][0] =  sinxp * sinyp;
            pm[1][1] =  cosyp;
            pm[1][2] =  cosxp * sinyp;
            pm[2][0] =  sinxp * cosyp;
            pm[2][1] = -sinyp;
            pm[2][2] =  cosxp * cosyp;
            
            return pm;
    }
    
    public double [][] precession (double ttt) {
        double ttt2 = ttt * ttt;
        double ttt3 = ttt * ttt2;
        double[][] prec = new double[3][3];
        double zeta =             2306.2181*ttt + 0.30188*ttt2 + 0.017998*ttt3;
        zeta = arcsec2rad(zeta);
        double theta=             2004.3109*ttt - 0.42665*ttt2 - 0.041833*ttt3;
        theta = arcsec2rad(theta);
        double z    =             2306.2181*ttt + 1.09468*ttt2 + 0.018203*ttt3;
        z = arcsec2rad(z);
        
        double coszeta  = Math.cos(zeta);
        double sinzeta  = Math.sin(zeta);
        double costheta = Math.cos(theta);
        double sintheta = Math.sin(theta);
        double cosz     = Math.cos(z);
        double sinz     = Math.sin(z);
        
            prec[0][0] =  coszeta * costheta * cosz - sinzeta * sinz;
            prec[0][1] =  coszeta * costheta * sinz + sinzeta * cosz;
            prec[0][2] =  coszeta * sintheta;
            prec[1][0] = -sinzeta * costheta * cosz - coszeta * sinz;
            prec[1][1] = -sinzeta * costheta * sinz + coszeta * cosz;
            prec[1][2] = -sinzeta * sintheta;
            prec[2][0] = -sintheta * cosz;
            prec[2][1] = -sintheta * sinz;
            prec[2][2] =  costheta;
        return prec;
    }
    
    public double[][] nutation(double ttt, double ddpsi, double ddeps) {
        double[][] nut = new double[3][3];
        ddpsi = arcsec2rad(ddpsi);
        ddeps = arcsec2rad(ddeps);
        double ttt2 = ttt * ttt;
        double ttt3 = ttt * ttt2;
        meaneps = -46.8150 *ttt - 0.00059 *ttt2 + 0.001813 *ttt3 + 84381.448;
        meaneps = arcsec2rad(meaneps);
        meaneps = remainder(meaneps, 2*PI);
        
        double l = ( (((0.064) * ttt + 31.310) * ttt + 1717915922.6330 ) * ttt) / 3600.0 + 134.96298139;
        l = remainder(l, 360.0);
        l = Math.toRadians(l);
    	double 	l1 = ( (((-0.012) * ttt - 0.577) * ttt + 129596581.2240 ) * ttt) / 3600.0 + 357.52772333;
        l1 = remainder(l1, 360.0);
        l1 = Math.toRadians(l1);
        double 	f = ( (((0.011) * ttt - 13.257) * ttt + 1739527263.1370 ) * ttt) / 3600.0 + 93.27191028;
        f = remainder(f, 360.0);
        f = Math.toRadians(f);
        double  d = ( (((0.019) * ttt - 6.891) * ttt + 1602961601.3280 ) * ttt) / 3600.0 + 297.85036306;
        d = remainder(d, 360.0);
        d = Math.toRadians(d);
        omega = ( (((0.008) * ttt + 7.455) * ttt - 6962890.5390 ) * ttt) / 3600.0 + 125.04452222;
        omega = remainder(omega, 360.0);
        omega = Math.toRadians(omega);
        
        double trueeps = 0; // true obliquity of the ecliptic
        double deltaeps = 0; // change in obliquity
        double tempval;
        for (int i = 105; i >= 0; i--) {
            tempval= iar80[i][0]*l + iar80[i][1]*l1 + iar80[i][2]*f + iar80[i][3]*d + iar80[i][4]*omega;
            deltapsi= deltapsi + (rar80[i][0]+rar80[i][1]*ttt) * Math.sin( tempval );
            deltaeps= deltaeps + (rar80[i][2]+rar80[i][3]*ttt) * Math.cos( tempval ); 
        }
        
        deltapsi = remainder((deltapsi + ddpsi), 2*PI);
        deltaeps = remainder((deltaeps + ddeps), 2*PI);
        trueeps  = meaneps + deltaeps;
        
        double cospsi  = Math.cos(deltapsi);
        double sinpsi  = Math.sin(deltapsi);
        double coseps  = Math.cos(meaneps);
        double sineps  = Math.sin(meaneps);
        double costrueeps = Math.cos(trueeps);
        double sintrueeps = Math.sin(trueeps);
        
        nut[0][0] =  cospsi;
        nut[0][1] =  costrueeps * sinpsi;
        nut[0][2] =  sintrueeps * sinpsi;
        nut[1][0] = -coseps * sinpsi;
        nut[1][1] =  costrueeps * coseps * cospsi + sintrueeps * sineps;
        nut[1][2] =  sintrueeps * coseps * cospsi - sineps * costrueeps;
        nut[2][0] = -sineps * sinpsi;
        nut[2][1] =  costrueeps * sineps * cospsi - sintrueeps * coseps;
        nut[2][2] =  sintrueeps * sineps * cospsi + costrueeps * coseps;
        
        return nut;
    }
    
    public static double remainder(double a, double b) {
        // this method returns the indentical remainder as of Matlab function rem(a, b), with
        // description as follows:
        // rem(a,b) returns the remainder after division of a by b, where a is the dividend and b is the divisor. 
        // This function is often called the remainder operation, which can be expressed as r = a - b.*fix(a./b). 
        // The rem function follows the convention that rem(a,0) is NaN.
                
        return a - b * roundToZero(a,b);
    }
    
    public static double roundToZero (double a, double b) {
        if (a/b < 0) {
            return Math.ceil(a/b);
        }
        return Math.floor(a/b);
    }
    
    public double[][] sidereal(double jdut1) {
        double[][] st = new double[3][3];
        double tut1 = ( jdut1 - 2451545.0 ) / 36525.0; 
        // Greenwich Mean Sidereal Time
        double GMST = - 6.2e-6 * tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 + (876600.0 * 3600.0 + 8640184.812866) * tut1 + 67310.54841;
        GMST = Math.toRadians(GMST); // rad
        GMST = GMST / 240.0;  //(360.0 / 86400.0); // GST in deg
        GMST = remainder(GMST, 2*PI);
        if (GMST < 0) { // ensure GST is within [0, 2*pi]
            GMST += 2 * PI;
        }
//        GMST = 5.459562587690165;
//        System.out.println("    GMST: " + Math.toDegrees(GMST));
        double AST; // Apprent Sidereal Time
        // after 1997, kinematic terms apply as well as gemoetric in eqe
        if ((jdut1 > 2450449.5 ) && (eqeterms > 0)) {
            AST = GMST + deltapsi* Math.cos(meaneps) 
                + arcsec2rad(0.00264) * Math.sin(omega) 
                + arcsec2rad(0.000063) * Math.sin(2.0*omega);
        }
        else {
            AST= GMST + deltapsi* Math.cos(meaneps);
        }
        AST = remainder(AST, 2*PI);
//        System.out.println("    AST: " + Math.toDegrees(AST));

        double thetasa = 7.29211514670698e-05 * (1.0  - lod/86400.0 );
        double omegaEarth = thetasa;
        
        st = GetC3(AST);
        
        // sidereal time rate matrix
        stdot[0][0] = -omegaEarth * Math.sin(AST);
        stdot[0][1] = -omegaEarth * Math.cos(AST);
        stdot[0][2] =  0.0;
        stdot[1][0] =  omegaEarth * Math.cos(AST);
        stdot[1][1] = -omegaEarth * Math.sin(AST);
        stdot[1][2] =  0.0;
        stdot[2][0] =  0.0;
        stdot[2][1] =  0.0;
        stdot[2][2] =  0.0;
        
        return st;
    }
    
    public static double arcsec2rad(double a) {
        return Math.toRadians(a/3600);
    }
    
    public static double[][] GetC1(double angle) {
        double[][] C1 = new double[3][3];
        double c = Math.cos(angle);
        double s = Math.sin(angle);

        C1[0][0] = 1;
        C1[0][1] = 0;
        C1[0][2] = 0;
        C1[1][0] = 0;
        C1[1][1] = c;
        C1[1][2] = -s;
        C1[2][0] = 0;
        C1[2][1] = s;
        C1[2][2] = c;
        return C1;
    }
    
    public static double[][] GetC2(double angle) {
        double[][] C2 = new double[3][3];
        double c = Math.cos(angle);
        double s = Math.sin(angle);

        C2[0][0] = c;
        C2[0][1] = 0;
        C2[0][2] = s;
        C2[1][0] = 0;
        C2[1][1] = 1;
        C2[1][2] = 0;
        C2[2][0] = -s;
        C2[2][1] = 0;
        C2[2][2] = c;
        return C2;
    }
    
    public static double[][] GetC3(double angle) {
        double[][] C3 = new double[3][3];
        double c = Math.cos(angle);
        double s = Math.sin(angle);

        C3[0][0] = c;
        C3[0][1] = -s;
        C3[0][2] = 0;
        C3[1][0] = s;
        C3[1][1] = c;
        C3[1][2] = 0;
        C3[2][0] = 0;
        C3[2][1] = 0;
        C3[2][2] = 1;
        return C3;
    }
    
    public static double[] VectorAdd(double[] a, double[] b) {
        int aRows = a.length;
        int bRows = b.length;
        double[] c = new double[aRows];
        if (aRows == bRows) {
            for (int k = 0; k < aRows; k++) {
                c[k] = a[k] + b[k];
            }
        } else {
            System.out.println("Array addition where the number of A elements (" + aRows + ") does not match the number of B elements (" + bRows + ")");
        }
        return c;
    }

    public static double[] matrixMult(double[][] a, double[] b) {
        int aRows = a.length;
        int aCols = a[0].length;
        int bRows = b.length;
        double[] c = new double[aRows];
        if (aCols == bRows) {
            for (int i = 0; i < aRows; i++) {
                for (int k = 0; k < aCols; k++) {
                    c[i] += a[i][k] * b[k];
                }
            }
        } else {
            System.out.println("Matrix multiplication where the number of aCols (" + aCols + ") does not match the number of bRows (" + bRows + ")");
        }
        return c;
    }

    public static double[][] matrixMult(double[][] a, double[][] b) {
        int aRows = a.length;
        int aCols = a[0].length;
        int bRows = b.length;
        int bCols = b[0].length;
        double[][] c = new double[aRows][bCols];
        if (aCols == bRows) {
            for (int i = 0; i < aRows; i++) {
                for (int j = 0; j < bCols; j++) {
                    c[i][j] = 0.0;
                    for (int k = 0; k < aCols; k++) {
                        c[i][j] += a[i][k] * b[k][j];
                    }
                }
            }
        } else {
            System.out.println("Maxtrix multiplication where the number of aCols " + aCols + ") does not match the number of bRows (" + bRows
                    + ")");
        }
        return c;
    }

    public static double[][] matrixMult(double[][] a, double b) {
        int aRows = a.length;
        int aCols = a[0].length;
        double[][] ret = new double[aRows][aCols];
        for (int r = 0; r < aRows; r++) {
            for (int c = 0; c < aCols; c++) {
                ret[r][c] = a[r][c] * b;
            }
        }
        return ret;
    }
    
    public static double[] crossProduct(double[] a, double[] b) {
		double[] c = {0, 0, 0};
		c[0] = a[1] * b[2] - a[2] * b[1];
		c[1] = a[2] * b[0] - a[0] * b[2];
		c[2] = a[0] * b[1] - a[1] * b[0];
		return c;
    }

    public static double[][] transpose(double[][] orig) {
        double[][] ret = new double[orig[0].length][orig.length];
        for (int rows = 0; rows < orig.length; rows++) {
            for (int cols = 0; cols < orig[0].length; cols++) {
                ret[cols][rows] = orig[rows][cols];
            }
        }
        return ret;
    }

    public static List<String> getFileContent(File f) throws IOException {
        List<String> ls = FileUtils.readLines(f, "utf-8");
        return ls;
    }
    
    public static void main(String[] args) {
        double a = -3.0;
        double b = -11.0;
        System.out.println(remainder(a, b));
        System.out.println(roundToZero(a, b));
        System.out.println(a % b);

//        double[] r_ecef = {-5667.275927,  1909.624625,  -3932.869963};
//        double[] r_ecef = {2517.147128,  4906.513829,  -4562.952512};
//        double[] r_ecef = {-1033.4793830,  7901.2952754,  6380.3565958};
        double[] r_ecef = {4706.438724, -4290.419904, 3274.510499};
//        double[] v_ecef = {3.764433,  3.317829,   5.638725};
//        double[] v_ecef = {-3.225636520,  -2.872451450,   5.531924446};
        double[] v_ecef = {1.341008, -3.530559, -6.537492};

        try {
            ecef2ecij2000 eci = new ecef2ecij2000(r_ecef, v_ecef, "2018-285-22:15:00.055");
//            System.out.println(Arrays.toString(eci.r_eci));
//            System.out.println("");
//            System.out.println(Arrays.toString(eci.v_eci));
        } catch (ParseException ex) {
            Logger.getLogger(ecef2ecij2000.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
