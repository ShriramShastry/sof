#!/bin/bash
# SPDX-License-Identifier: BSD-3-Clause
# Copyright(c) 2021 Google Inc. All rights reserved.

cmake -DCMAKE_INSTALL_PREFIX=install
make install -j $(nproc)
