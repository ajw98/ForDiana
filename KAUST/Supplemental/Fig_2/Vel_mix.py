#!/usr/bin/env python3

from PIL import Image

images = [
    Image.open("11/plt/nfields_00026.png"),
    Image.open("11/plt/nfields_00046.png"),
    Image.open("12/plt/nfields_00050.png"),
    Image.open("12/plt/nfields_00068.png"),
    Image.open("11/vel/nvel_00026.png"),
    Image.open("11/vel/nvel_00046.png"),
    Image.open("12/vel/nvel_00050.png"),
    Image.open("12/vel/nvel_00068.png")
    ]


for i in range(4, 8):
    img = images[i]
    images[i] = img.crop((
        0,              # left
        200,            # top
        img.width,      # right
        img.height      # bottom
    ))


# Split into two rows
row1 = images[:4]
row2 = images[4:]

# Width of each column
col_widths = [
    max(row1[0].width, row2[0].width),
    max(row1[1].width, row2[1].width),
    max(row1[2].width, row2[2].width),
    max(row1[3].width, row2[3].width),
    ]

# Height of each row
row1_height = max(img.height for img in row1)
row2_height = max(img.height for img in row2)

total_width = sum(col_widths)
total_height = row1_height + row2_height
print(total_width, total_height)
combined = Image.new(
    "RGB",
    (total_width, total_height),
    "white"
)

# First row
x = 0
for i, img in enumerate(row1):
    combined.paste(img, (x, 0))
    x += col_widths[i]

# Second row
x = 0
for i, img in enumerate(row2):
    combined.paste(img, (x, row1_height))
    x += col_widths[i]

combined.save("Combined2.png", dpi=(600, 600))
combined.save("Combined2.pdf", dpi=(600, 600))
