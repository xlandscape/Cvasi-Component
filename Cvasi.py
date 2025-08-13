import numpy
import os
import h5py
import datetime
import base
import attrib


class Cvasi(base.Component):
    def __init__(self, name, observer, store):
        super(Cvasi, self).__init__(name, observer, store)
        self._inputs = base.InputContainer(self, [
            base.Input(
                "ProcessingPath",
                (attrib.Class(str), attrib.Unit(None), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "k_phot_fix",
                (attrib.Class(bool), attrib.Unit(None), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "k_phot_max",
                (attrib.Class(float), attrib.Unit("1/d"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "k_resp",
                (attrib.Class(float), attrib.Unit("1/d"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "k_loss",
                (attrib.Class(float), attrib.Unit("1/d"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "mass_per_frond",
                (attrib.Class(float), attrib.Unit("g"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "BMw2BMd",
                (attrib.Class(float), attrib.Unit("1"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "Emax",
                (attrib.Class(float), attrib.Unit("1"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "EC50",
                (attrib.Class(float), attrib.Unit("µg/L"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "b",
                (attrib.Class(float), attrib.Unit("1"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "P_up",
                (attrib.Class(float), attrib.Unit("cm/d"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "AperBM",
                (attrib.Class(float), attrib.Unit("cm²/g"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "Kbm",
                (attrib.Class(float), attrib.Unit("1"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "P_Temp",
                (attrib.Class(bool), attrib.Unit(None), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "MolWeight",
                (attrib.Class(float), attrib.Unit("g/mol"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "Tmin",
                (attrib.Class(float), attrib.Unit("°C"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "Tmax",
                (attrib.Class(float), attrib.Unit("°C"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "Topt",
                (attrib.Class(float), attrib.Unit("°C"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "t_ref",
                (attrib.Class(float), attrib.Unit("°C"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "Q10",
                (attrib.Class(float), attrib.Unit("1"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "k_0",
                (attrib.Class(float), attrib.Unit("1/d"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "a_k",
                (attrib.Class(float), attrib.Unit("m²*d/kJ"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "C_P",
                (attrib.Class(float), attrib.Unit("mg/L"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "CP50",
                (attrib.Class(float), attrib.Unit("mg/L"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "a_P",
                (attrib.Class(float), attrib.Unit("1"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "KiP",
                (attrib.Class(float), attrib.Unit("mg/L"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "C_N",
                (attrib.Class(float), attrib.Unit("mg/L"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "CN50",
                (attrib.Class(float), attrib.Unit("mg/L"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "a_N",
                (attrib.Class(float), attrib.Unit("1"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "KiN",
                (attrib.Class(float), attrib.Unit("mg/L"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "BM50",
                (attrib.Class(float), attrib.Unit("g/m²"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "BM",
                (attrib.Class(float), attrib.Unit("g/m²"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "E",
                (attrib.Class(float), attrib.Unit("1"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "M_int",
                (attrib.Class(float), attrib.Unit("µg/m²"), attrib.Scales("global")),
                self.default_observer
            ),
            base.Input(
                "Temperature",
                (attrib.Class(numpy.ndarray), attrib.Unit("°C"), attrib.Scales("time/day")),
                self.default_observer
            ),
            base.Input(
                "Radiation",
                (attrib.Class(numpy.ndarray), attrib.Unit("kJ/(m²*d)"), attrib.Scales("time/day")),
                self.default_observer
            ),
            base.Input(
                "Concentration",
                (attrib.Class(numpy.ndarray), attrib.Unit("mg/m³"), attrib.Scales("time/hour, space/reach")),
                self.default_observer
            )
        ])
        self._outputs = base.OutputContainer(
            self,
            [
                base.Output("BM_out", store, self),
                base.Output("E_out", store, self),
                base.Output("M_int_out", store, self),
                base.Output("C_int_out", store, self),
                base.Output("FrondNo_out", store, self),
                base.Output("BM.EP10", store, self),
                base.Output("r.EP10", store, self),
                base.Output("BM.EP50", store, self),
                base.Output("r.EP50", store, self)
            ]
        )

    def run(self) -> None:
        processing_path = self.inputs["ProcessingPath"].read().values
        os.makedirs(processing_path)
        concentrations_description = self.inputs["Concentration"].describe()
        number_days = int(concentrations_description["shape"][0] / 24)
        h5file = os.path.join(processing_path, "cvasi.h5")
        with h5py.File(h5file, "a") as f:
            f.create_dataset("/init/BM", (1,), "d", self.inputs["BM"].read().values)
            f.create_dataset("/init/M_int", (1,), "d", self.inputs["M_int"].read().values)
            f.create_dataset("/init/E", (1,), "d", self.inputs["E"].read().values)
            f.create_dataset("/times", (1,), "i", number_days)
            f.create_dataset("/param/k_phot_fix", (1,), "i", 1 if self.inputs["k_phot_fix"].read().values else 0)
            f.create_dataset("/param/k_phot_max", (1,), "d", self.inputs["k_phot_max"].read().values)
            f.create_dataset("/param/k_resp", (1,), "d", self.inputs["k_resp"].read().values)
            f.create_dataset("/param/k_loss", (1,), "d", self.inputs["k_loss"].read().values)
            f.create_dataset("/param/mass_per_frond", (1,), "d", self.inputs["mass_per_frond"].read().values)
            f.create_dataset("/param/BMw2BMd", (1,), "d", self.inputs["BMw2BMd"].read().values)
            f.create_dataset("/param/Emax", (1,), "d", self.inputs["Emax"].read().values)
            f.create_dataset("/param/EC50", (1,), "d", self.inputs["EC50"].read().values)
            f.create_dataset("/param/b", (1,), "d", self.inputs["b"].read().values)
            f.create_dataset("/param/P_up", (1,), "d", self.inputs["P_up"].read().values)
            f.create_dataset("/param/AperBM", (1,), "d", self.inputs["AperBM"].read().values)
            f.create_dataset("/param/Kbm", (1,), "d", self.inputs["Kbm"].read().values)
            f.create_dataset("/param/P_Temp", (1,), "d", self.inputs["P_Temp"].read().values)
            f.create_dataset("/param/MolWeight", (1,), "d", self.inputs["MolWeight"].read().values)
            f.create_dataset("/param/Tmin", (1,), "d", self.inputs["Tmin"].read().values)
            f.create_dataset("/param/Tmax", (1,), "d", self.inputs["Tmax"].read().values)
            f.create_dataset("/param/Topt", (1,), "d", self.inputs["Topt"].read().values)
            f.create_dataset("/param/t_ref", (1,), "d", self.inputs["t_ref"].read().values)
            f.create_dataset("/param/Q10", (1,), "d", self.inputs["Q10"].read().values)
            f.create_dataset("/param/k_0", (1,), "d", self.inputs["k_0"].read().values)
            f.create_dataset("/param/a_k", (1,), "d", self.inputs["a_k"].read().values)
            f.create_dataset("/param/C_P", (1,), "d", self.inputs["C_P"].read().values)
            f.create_dataset("/param/CP50", (1,), "d", self.inputs["CP50"].read().values)
            f.create_dataset("/param/a_P", (1,), "d", self.inputs["a_P"].read().values)
            f.create_dataset("/param/KiP", (1,), "d", self.inputs["KiP"].read().values)
            f.create_dataset("/param/C_N", (1,), "d", self.inputs["C_N"].read().values)
            f.create_dataset("/param/CN50", (1,), "d", self.inputs["CN50"].read().values)
            f.create_dataset("/param/a_N", (1,), "d", self.inputs["a_N"].read().values)
            f.create_dataset("/param/KiN", (1,), "d", self.inputs["KiN"].read().values)
            f.create_dataset("/param/BM50", (1,), "d", self.inputs["BM50"].read().values)
            f.create_dataset(
                "/envir/tmp",
                dtype="d",
                data=self.inputs["Temperature"].read(
                    select={
                        "time/day": {
                            "from": concentrations_description["offsets"][0].date(),
                            "to": (
                                    concentrations_description["offsets"][0] +
                                    datetime.timedelta(days=1, hours=concentrations_description["shape"][0])
                            ).date()
                        }
                    }
                ).values
            )
            f.create_dataset(
                "/envir/irr",
                dtype="d",
                data=self.inputs["Radiation"].read(
                    select={
                        "time/day": {
                            "from": concentrations_description["offsets"][0].date(),
                            "to": (
                                    concentrations_description["offsets"][0] +
                                    datetime.timedelta(days=1, hours=concentrations_description["shape"][0])
                            ).date()
                        }
                    }
                ).values
            )
            concentrations_dataset = f.create_dataset(
                "envir/conc",
                concentrations_description["shape"],
                "d",
                compression="gzip",
                chunks=(concentrations_description["shape"][0], 1)
            )
            for reach in range(concentrations_description["shape"][1]):
                concentrations_dataset[:, reach:reach + 1] = self.inputs["Concentration"].read(
                    slices=(slice(concentrations_description["shape"][0]), slice(reach, reach + 1))).values
            for output in (("BM", "g/m²"), ("E", "1"), ("M_int", "µg/m²"), ("C_int", "µg/L"), ("FrondNo", 1)):
                f.create_dataset(
                    f"out/{output[0]}",
                    (number_days, concentrations_description["shape"][1]),
                    "d",
                    compression="gzip",
                    chunks=(number_days, 1)
                )
                self.outputs[f"{output[0]}_out"].set_values(
                    numpy.ndarray,
                    shape=(number_days, concentrations_description["shape"][1]),
                    data_type="d",
                    chunks=(number_days, 1),
                    scales="time/day, space/reach",
                    unit=output[1],
                    element_names=(None, concentrations_description["element_names"][1]),
                    offset=(concentrations_description["offsets"][0].date(), None),
                    geometries=(None, concentrations_description["geometries"][1])
                )
            for output in ("BM.EP10", "r.EP10", "BM.EP50", "r.EP50"):
                f.create_dataset(f"effects/{output}", (concentrations_description["shape"][1],), "d")
        r_exe = os.path.join(os.path.dirname(__file__), "module", "R-4.5.1", "bin", "x64", "Rscript.exe")
        r_script = os.path.join(os.path.dirname(__file__), "module", "cvasi.R")
        library_path = os.path.join(os.path.dirname(__file__), "module", "R-4.3.2", "library")
        base.run_process(
            (r_exe, "--vanilla", r_script, h5file),
            processing_path,
            self.default_observer,
            {"R_LIBS": library_path, "R_LIBS_USER": library_path}
        )
        with h5py.File(h5file) as f:
            for output in ("BM", "E", "M_int", "C_int", "FrondNo"):
                for reach in range(concentrations_description["shape"][1]):
                    self.outputs[f"{output}_out"].set_values(
                        f[f"out/{output}"][:, reach:(reach + 1)],
                        slices=(slice(number_days), slice(reach, reach + 1)),
                        create=False
                    )
            for output in ("BM.EP10", "r.EP10", "BM.EP50", "r.EP50"):
                self.outputs[f"{output}"].set_values(
                    f[f"effects/{output}"][:],
                    shape=(concentrations_description["shape"][1],),
                    data_type="d",
                    chunks=(concentrations_description["shape"][1],),
                    scales="space/reach",
                    unit="%",
                    element_names=(concentrations_description["element_names"][1],),
                    offset=(None,),
                    geometries=(concentrations_description["geometries"][1],)
                )
