"""Motor and generator energy and sign checks for the RM3 hydraulic PTO."""

import numpy as np
import pytest

from wecsim.electricGenerator import EquivalentCircuitGenerator
from wecsim.hydraulic import ConstantEfficiencyHydraulicMotor


def test_motor_power_ratio_matches_both_efficiencies():
    motor = ConstantEfficiencyHydraulicMotor(120e-6, .9, .85)
    pressure_drop = 12e6
    speed = 100.
    hydraulic_power = pressure_drop * motor.flow(speed)
    shaft_power = motor.torque(pressure_drop) * speed
    np.testing.assert_allclose(shaft_power / hydraulic_power, .9 * .85)


def test_generator_drive_and_back_emf_signs():
    generator = EquivalentCircuitGenerator(.8, .8, .8, .8, .8)
    assert generator.speed_rate(0, 0, 160) == 200
    assert generator.current_rate(100, 0, 0) == -100
    assert generator.electromagnetic_torque(-10) == -8


def test_drive_rejects_invalid_efficiency_and_inertia():
    with pytest.raises(ValueError, match="efficiencies"):
        ConstantEfficiencyHydraulicMotor(120e-6, 1.1, .85)
    with pytest.raises(ValueError, match="invalid"):
        EquivalentCircuitGenerator(.8, .8, .8, 0, .8)
