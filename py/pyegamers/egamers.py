import numpy as np
import matplotlib.pyplot as plt


class EGAMERS:

    from .ncfield import (
        read_ncfield,
        compute_ncfield,
        compute_gamma,
        plot_field_radial,
        plot_field_time,
        plot_spectrogram,
        plot_field_energy
    )

    class OutputData:
        pass

    def __init__(self, folder_path="./"):
        """! Read data from given path"""
        import os

        abspath = os.path.abspath(folder_path)

        namelist_name = os.path.join(abspath, "namelist.in")
        if os.path.exists(namelist_name):
            import f90nml
            self.namelist = f90nml.read(namelist_name)

        ncfield_name = os.path.join(abspath, "snapshot.field.nc")
        if os.path.exists(ncfield_name):
            self.ncfield = self.read_ncfield(ncfield_name)

        fqc_name = os.path.join(abspath, "fqc.out")
        if os.path.exists(fqc_name):
            self.fqcdata = self.read_fqc(fqc_name)

        field_name = os.path.join(abspath, "field.out")
        if os.path.exists(field_name):
            self.eigenfield = self.read_field(field_name)

        map_name = os.path.join(abspath, "map.out")
        if os.path.exists(map_name):
            self.fqcmapdata = self.read_fqcmap(map_name)

        orbit_name = os.path.join(abspath, "orbit.out")
        if os.path.exists(orbit_name):
            self.orbitdata = self.read_orbit(orbit_name)

    def read_fqc(self, fqc_name):
        fo = self.OutputData()

        f = open(fqc_name)
        line1 = f.readline()
        ns = list(map(int, line1.split()))
        fo.neigen = ns[0]
        fo.nr = ns[1]
        fo.nregam = ns[2]

        fo.r = np.zeros((fo.nr), dtype=float)
        fo.omglocal = np.zeros((fo.nr), dtype=float)
        fo.omgegamlocal = np.zeros((fo.nregam), dtype=float)
        fo.regamlocal = np.zeros((fo.nregam), dtype=float)
        fo.gammaegamlocal = np.zeros((fo.nregam), dtype=float)
        fo.omgglobal = np.zeros((fo.neigen), dtype=float)
        fo.gammaglobal = np.zeros((fo.neigen), dtype=float)

        for i1 in range(fo.neigen):
            line = f.readline()
            number_r = list(map(float, line.split()))
            fo.omgglobal[i1] = number_r[0]
            fo.gammaglobal[i1] = number_r[1]

        for i2 in range(fo.nr):
            line = f.readline()
            number_r = list(map(float, line.split()))
            fo.r[i2] = number_r[0]
            fo.omglocal[i2] = number_r[1]

        for i3 in range(fo.nregam):
            line = f.readline()
            number_r = list(map(float, line.split()))
            fo.regamlocal[i3] = number_r[0]
            fo.omgegamlocal[i3] = number_r[1]
            fo.gammaegamlocal[i3] = number_r[2]

        return fo

    def read_field(self, field_name):
        fo = self.OutputData()

        f = open(field_name)
        line1 = f.readline()
        ns = list(map(int, line1.split()))
        fo.neigen = ns[0]
        fo.nr = ns[1]

        fo.r = np.zeros((fo.neigen, fo.nr), dtype=float)
        fo.err = np.zeros((fo.neigen, fo.nr), dtype=float)
        fo.eri = np.zeros((fo.neigen, fo.nr), dtype=float)

        for i1 in range(fo.neigen):
            for i2 in range(fo.nr):
                line = f.readline()
                number_r = list(map(float, line.split()))
                fo.r[i1, i2] = number_r[0]
                fo.err[i1, i2] = number_r[1]
                fo.eri[i1, i2] = number_r[2]

        return fo

    def read_fqcmap(self, fqcmap_name=None):
        mo = self.OutputData()

        f = open(fqcmap_name)
        line1 = f.readline()
        headerline = list(map(int, line1.split()))
        mo.itrap = headerline[0]
        mo.icop = headerline[1]
        mo.ictp = headerline[2]

        if mo.itrap == 1:
            line2 = f.readline()
            trapline = list(map(int, line2.split()))
            mo.ipphin = trapline[0]
            mo.neen = trapline[1]

            mo.pphigridn = np.zeros((mo.ipphin, mo.neen), dtype=float)
            mo.een = np.zeros((mo.ipphin, mo.neen), dtype=float)
            mo.fqcn = np.zeros((mo.ipphin, mo.neen), dtype=float)

            for i1 in range(mo.ipphin):
                for i2 in range(mo.neen):
                    line = f.readline()
                    data1 = list(map(float, line.split()))
                    mo.pphigridn[i1, i2] = data1[0]
                    mo.een[i1, i2] = data1[1]
                    mo.fqcn[i1, i2] = data1[2]

        return mo

    def read_orbit(self, orbit_name):
        oo = self.OutputData()

        f = open(orbit_name)
        line1 = f.readline()
        omgline = list(map(float, line1.split()))
        omega = omgline[0]
        line2 = f.readline()
        ns = list(map(int, line2.split()))
        n = ns[0]

        r = np.array((0.0,) * n)
        theta = np.array((0.0,) * n)

        for i1 in range(n):
            line = f.readline()
            data = list(map(float, line.split()))
            r[i1] = data[0]
            theta[i1] = data[1]

        x = r * np.cos(theta)
        y = r * np.sin(theta)

        thetam = np.arange(0.0, 2 * np.pi, 2 * np.pi / n)
        xm = np.cos(thetam)
        ym = np.sin(thetam)

        oo.x = x
        oo.y = y
        oo.thetam = thetam
        oo.xm = xm
        oo.ym = ym
        oo.omega = omega
        oo.n = n

        return oo

    def plot_fqcmap(self, fqcmapdata=None):
        if fqcmapdata == None:
            fqcmapdata = self.fqcmapdata

            plt.figure(1)
            CS = plt.contourf(
                fqcmapdata.pphigridn.reshape(fqcmapdata.ipphin, fqcmapdata.neen),
                fqcmapdata.een.reshape(fqcmapdata.ipphin, fqcmapdata.neen),
                fqcmapdata.fqcn.reshape(fqcmapdata.ipphin, fqcmapdata.neen),
                cmap=plt.cm.bone,
            )
            cbar = plt.colorbar(CS)
            cbar.ax.set_ylabel(r"$\omega_b$")
            plt.title("Trapped particle frequency map")
            plt.xlabel(r"$P_\varphi / e \Psi$")
            plt.ylabel("E (keV)")
            plt.show()

    def plot_fqc(self, fqcdata=None):
        if fqcdata == None:
            fqcdata = self.fqcdata

        (ll,) = plt.plot(fqcdata.r, fqcdata.omglocal)
        (llegam,) = plt.plot(fqcdata.regamlocal, fqcdata.omgegamlocal)
        for i1 in range(fqcdata.neigen):
            (gg,) = plt.plot([0, 1], [fqcdata.omgglobal[i1], fqcdata.omgglobal[i1]])
        #    leg = plt.legend([gg], ['Global Mode ' + str(i1+1)])
        #    ax = plt.gca().add_artist(leg)

        plt.legend([ll, llegam, gg], ["GAM continuum", "EGAM continuum", "Global EGAM"])
        plt.ylabel(r"$Re(\Omega)$")
        plt.xlabel("r")
        plt.title(
            r"$\gamma/\omega$ = " + str(fqcdata.gammaglobal[0] / fqcdata.omgglobal[0])
        )

        plt.figure(2)
        (llegam,) = plt.plot(fqcdata.regamlocal, fqcdata.gammaegamlocal)
        for i1 in range(fqcdata.neigen):
            (gg,) = plt.plot([0, 1], [fqcdata.gammaglobal[i1], fqcdata.gammaglobal[i1]])
        plt.legend([llegam, gg], ["EGAM continuum", "Global EGAM"])
        plt.ylabel(r"$\gamma=Im(\Omega)$")
        plt.xlabel("r")
        plt.title(
            r"$\gamma/\omega$ = " + str(fqcdata.gammaglobal[0] / fqcdata.omgglobal[0])
        )
        plt.show()

    def plot_field(self, fielddata=None):
        if fielddata == None:
            fielddata = self.eigenfield

        for i1 in range(fielddata.neigen):
            r = fielddata.r[i1]
            err = fielddata.err[i1]
            eri = fielddata.eri[i1]

            era = np.sqrt(np.square(err) + np.square(eri))
            plt.figure(i1)
            (rep,) = plt.plot(r, err)
            (imp,) = plt.plot(r, eri)
            (abp,) = plt.plot(r, era)
            plt.ylabel("Er")
            plt.xlabel("r")
            plt.legend(
                [rep, imp, abp],
                [
                    "Real (Mode " + str(i1 + 1) + ")",
                    "Imag (Mode " + str(i1 + 1) + ")",
                    "Abs (Mode " + str(i1 + 1) + ")",
                ],
            )
            plt.show()

    def plot_orbit(self, orbitdata=None):
        if orbitdata == None:
            orbitdata = self.orbitdata

        plt.figure(1)
        plt.plot(orbitdata.x, orbitdata.y, orbitdata.xm, orbitdata.ym)
        plt.axis("equal")
        plt.title(r"$\omega_b$ = " + str(orbitdata.omega))
        plt.show()
