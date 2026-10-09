"""Planar linkage used by the published WECCCOMP WaveStar model."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class WaveStarLinkage:
    """Rotate the float and arm about A and measure the B-to-C PTO stroke.

    Points are ``(x, z)`` world coordinates in the neutral pose. The PTO
    stroke is positive when the B-to-C distance shortens, as in the source.
    """

    pivot: tuple[float, float] = (-0.438, 0.302)
    pto_base: tuple[float, float] = (-0.438, 0.714)
    pto_arm: tuple[float, float] = (-0.6214398, 0.3816858)

    def point_position(self, neutral_point, angle):
        offset = np.asarray(neutral_point, dtype=float) - self.pivot
        cosine, sine = np.cos(angle), np.sin(angle)
        return np.asarray(self.pivot) + np.stack((
            cosine * offset[0] + sine * offset[1],
            -sine * offset[0] + cosine * offset[1],
        ), axis=-1)

    def point_velocity(self, neutral_point, angle, angular_speed):
        offset = np.asarray(neutral_point, dtype=float) - self.pivot
        cosine, sine = np.cos(angle), np.sin(angle)
        return np.asarray(angular_speed)[..., None] * np.stack((
            -sine * offset[0] + cosine * offset[1],
            -cosine * offset[0] - sine * offset[1],
        ), axis=-1)

    def pto_stroke(self, angle):
        initial_length = np.linalg.norm(
            np.asarray(self.pto_arm) - self.pto_base,
        )
        current_length = np.linalg.norm(
            self.point_position(self.pto_arm, angle) - self.pto_base,
            axis=-1,
        )
        return initial_length - current_length

    def pto_speed(self, angle, angular_speed):
        arm = self.point_position(self.pto_arm, angle) - self.pto_base
        arm_speed = self.point_velocity(self.pto_arm, angle, angular_speed)
        return -np.sum(arm * arm_speed, axis=-1) / np.linalg.norm(arm, axis=-1)
