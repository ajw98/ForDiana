#!/usr/bin/env python3

#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# ------------------
# Files
# ------------------

top_left_file = "pol_00005.png"
top_right_file = "pol_00006.png"

bottom_left_file = "vel_00005.png"
bottom_right_file = "vel_00006.png"

crop_top_pixels = 100

# ------------------
# Load images
# ------------------

img_tl = mpimg.imread(top_left_file)
img_tr = mpimg.imread(top_right_file)

img_bl = mpimg.imread(bottom_left_file)
img_br = mpimg.imread(bottom_right_file)

# Crop bottom row images
img_bl = img_bl[crop_top_pixels:, :, :]
img_br = img_br[crop_top_pixels:, :, :]

# ------------------
# Figure sizing
# ------------------

row1_ratio = max(
    img_tl.shape[0] / img_tl.shape[1],
    img_tr.shape[0] / img_tr.shape[1]
)

row2_ratio = max(
    img_bl.shape[0] / img_bl.shape[1],
    img_br.shape[0] / img_br.shape[1]
)

figure_width = 6
figure_height = figure_width * (row1_ratio + row2_ratio) / 2

# ------------------
# Create layout
# ------------------

fig, axes = plt.subplots(
    2, 2,
    figsize=(figure_width, figure_height),
    dpi=300,
    gridspec_kw={"hspace": 0, "wspace": 0}
)

images = [
    [img_tl, img_tr],
    [img_bl, img_br]
]

for r in range(2):
    for c in range(2):
        axes[r, c].imshow(images[r][c])
        axes[r, c].axis("off")

plt.subplots_adjust(
    left=0,
    right=1,
    top=1,
    bottom=0,
    hspace=0,
    wspace=0
)

plt.savefig(
    "combined_2x2.png",
    dpi=300,
    bbox_inches="tight",
    pad_inches=0
)

plt.close()
