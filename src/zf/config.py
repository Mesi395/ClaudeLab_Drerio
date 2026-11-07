"""Schemas for the YAML files in config/ and markers/. A file that fails here fails the tests."""

from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator

Condition = Literal["uninjured", "dpi_3", "dpi_7", "dpi_30"]
Role = Literal[
    "core",
    "dissociation_sensitive",
    "non_specific",
    "library_size_suspect",
    "ortholog_mapping_ours",
]


class Strict(BaseModel):
    model_config = ConfigDict(extra="forbid", populate_by_name=True)


class Library(Strict):
    library: str
    gsm: str
    geo_label: str
    condition: Condition
    chamber: Literal["whole", "atrium", "ventricle"]
    chemistry: Literal["10x_v2", "10x_v3"]
    study: str
    expected_cells: int


class Expected(Strict):
    cells: int
    genes: int
    libraries: int


class Dataset(Strict):
    schema_: Literal["zf.dataset/v1"] = Field(alias="schema")
    name: str
    citation: str
    annotation: str
    expected: Expected
    covariates: list[str]
    libraries: list[Library]


class MarkerGene(Strict):
    symbol: str
    ensembl: str = Field(pattern=r"^ENSDARG\d{11}$")
    role: Role
    source: str
    evidence: str = ""
    caveat: str = ""
    status: Literal["candidate", "accepted", "flagged", "rejected"] = "candidate"
    proposed_by: Literal["literature_review", "spatial_data", "agent", "expert"]
    flags: list[str] = []


class Markers(Strict):
    schema_: Literal["zf.markers/v1"] = Field(alias="schema")
    axis: str
    annotation: str
    genes: list[MarkerGene]


class Expectation(Strict):
    higher: Condition
    than: Condition
    status: Literal["proposed", "confirmed", "dropped"]
    basis: str


class Axis(Strict):
    name: str
    markers: str
    score_roles: list[Role]
    expect: list[Expectation] = []


class ControlGene(Strict):
    symbol: str
    ensembl: str


class Control(Strict):
    name: str
    purpose: str
    genes: list[ControlGene] = []
    from_roles: list[Role] = []


class Axes(Strict):
    schema_: Literal["zf.axes/v1"] = Field(alias="schema")
    axes: list[Axis]
    controls: list[Control]


class Check(Strict):
    id: str
    family: Literal["signal", "confound", "robustness", "report"]
    detects: str | None
    question: str
    applies_to: str | list[str]
    method: str
    flag_if: str
    reference: str | None = None


class Checks(Strict):
    schema_: Literal["zf.checks/v2"] = Field(alias="schema")
    status: Literal["proposed", "frozen"]
    frozen_at: str | None
    seeds: list[int] = Field(min_length=2)
    bootstrap: int
    covariates: dict[str, str]
    checks: list[Check]

    @model_validator(mode="after")
    def one_check_per_defect(self):
        detected = [c.detects for c in self.checks if c.family != "report"]
        if None in detected or len(detected) != len(set(detected)):
            raise ValueError(
                "every gating check must detect exactly one defect, and no two the same"
            )
        confounds = {c.id for c in self.checks if c.family == "confound" and c.id != "batch"}
        if confounds != set(self.covariates):
            raise ValueError("each confound check needs exactly one covariate, named by its id")
        return self


class Split(Strict):
    test: list[str]
    train: Literal["rest"]


class Hillclimb(Strict):
    schema_: Literal["zf.hillclimb/v1"] = Field(alias="schema")
    surface: str
    objective: str
    constraints: dict
    split: Split
    rounds: dict
    on_stall: str
    final_test: str


def _load(model, path):
    return model.model_validate(yaml.safe_load(Path(path).read_text()))


def load_dataset(path) -> Dataset:
    return _load(Dataset, path)


def load_markers(path) -> Markers:
    return _load(Markers, path)


def load_axes(path) -> Axes:
    return _load(Axes, path)


def load_checks(path) -> Checks:
    return _load(Checks, path)


def load_hillclimb(path) -> Hillclimb:
    return _load(Hillclimb, path)
