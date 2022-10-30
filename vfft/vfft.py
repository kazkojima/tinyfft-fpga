#!/usr/bin/env python3
#
# Copyright (c) 2022 Kaz Kojima <kkojima@rr.iij4u.or.jp>
# SPDX-License-Identifier: CERN-OHL-W-2.0

from amaranth import *
from amaranth import Signal, Module, Elaboratable, Memory
from amaranth.utils import log2_int
from amaranth.cli import main

from amlib.test import GatewareTestCase, sync_test_case

import numpy as np
from math import cos, sin, pi
from pprint import pformat

class FixedPointFFT(Elaboratable):
    def __init__(self,
                 bitwidth: int=18,
                 pts:      int=1024,
                 verbose:  bool=True) -> None:

        self.bitwidth = bitwidth
        self.pts      = pts
        self.stages   = log2_int(pts)

        self.in_i = Signal(signed(bitwidth))
        self.in_q = Signal(signed(bitwidth))
        self.out_real = Signal(signed(bitwidth))
        self.out_imag = Signal(signed(bitwidth))
        self.strobe_in = Signal()
        self.strobe_out = Signal()

        self.Wr = [int(cos(k*2*pi/pts)*(2**(bitwidth-1))) for k in range(pts)]
        self.Wi = [int(-sin(k*2*pi/pts)*(2**(bitwidth-1))) for k in range(pts)]

        # 4-term Blackman-Harris window function
        self.wF = [int((0.35875-0.48829*cos(k*2*pi/pts)+0.14128*cos(k*4*pi/pts)-0.01168*cos(k*6*pi/pts))*(2**(bitwidth-1))) for k in range(pts)]

        assert pts == 2**self.stages, f"Points {pts} must be 2**stages {self.stages}"

    def elaborate(self, platform) -> Module:
        m = Module()

        width = self.bitwidth
        bw = width + self.stages
        pts = self.pts
        xr = Memory(width=bw, depth=self.pts, name="xr")
        xi = Memory(width=bw, depth=self.pts, name="xi")
        yr = Memory(width=bw, depth=self.pts, name="yr")
        yi = Memory(width=bw, depth=self.pts, name="yi")

        wF = Array(Const(v, signed(width+1)) for v in self.wF)

        Wr = Array(Const(v, signed(width+1)) for v in self.Wr)
        Wi = Array(Const(v, signed(width+1)) for v in self.Wi)
        
        m.submodules.xr_rd = xr_rd = xr.read_port()
        m.submodules.xr_wr = xr_wr = xr.write_port()
        m.submodules.xi_rd = xi_rd = xi.read_port()
        m.submodules.xi_wr = xi_wr = xi.write_port()
        m.submodules.yr_rd = yr_rd = yr.read_port()
        m.submodules.yr_wr = yr_wr = yr.write_port()
        m.submodules.yi_rd = yi_rd = yi.read_port()
        m.submodules.yi_wr = yi_wr = yi.write_port()

        N = self.stages
        idx = Signal(N+1)
        revidx = Signal(N)
        m.d.comb += revidx.eq(Cat([idx.bit_select(i,1) for i in reversed(range(N))]))

        # Window
        wf = Signal(signed(width+1))
        i_cooked = Signal(signed(width))
        q_cooked = Signal(signed(width))
        m.d.comb += [
            wf.eq(wF[idx]),
            i_cooked.eq((wf*self.in_i) >> (width-1)),
            q_cooked.eq((wf*self.in_q) >> (width-1)),
        ]

        # FFT
        widx = Signal(N)
        stage = Signal(range(N+1))
        mask = Signal(signed(N))

        ar = Signal(signed(bw))
        ai = Signal(signed(bw))
        br = Signal(signed(bw))
        bi = Signal(signed(bw))

        # coefficients
        wr = Signal(signed(width+1))
        wi = Signal(signed(width+1))

        m.d.comb += widx.eq(idx & mask)
        m.d.comb += wr.eq(Wr[widx])
        m.d.comb += wi.eq(Wi[widx])

        # complex multiplication
        mrr = Signal(signed(bw))
        mii = Signal(signed(bw))
        mri = Signal(signed(bw))
        mir = Signal(signed(bw))
        bwr = Signal(signed(bw))
        bwi = Signal(signed(bw))

        m.d.comb += [
            mrr.eq((br * wr) >> (width-1)),
            mii.eq((bi * wi) >> (width-1)),
            mri.eq((br * wi) >> (width-1)),
            mir.eq((bi * wr) >> (width-1)),
            bwr.eq(mrr - mii),
            bwi.eq(mri + mir),
        ]

        # butterfly
        si = Signal(signed(bw))
        sr = Signal(signed(bw))
        di = Signal(signed(bw))
        dr = Signal(signed(bw))
        m.d.comb += [
            sr.eq(ar + bwr),
            si.eq(ai + bwi),
            dr.eq(ar - bwr),
            di.eq(ai - bwi),
        ]

        # Control FSM
        with m.FSM(reset="IDLE"):
            with m.State("IDLE"):
                m.d.sync += idx.eq(0)
                m.next = "WINDOW"

            with m.State("WINDOW"):
                m.d.sync += xr_wr.en.eq(0)
                m.d.sync += xi_wr.en.eq(0)
                with m.If(idx >= pts):
                    m.d.sync += stage.eq(0)
                    m.d.sync += idx.eq(0)
                    m.d.sync += mask.eq(~((2 << (N-2))-1))
                    m.next = "FFTLOOP"
                with m.Else():
                    m.next = "WINDOW_MUL"

            with m.State("WINDOW_MUL"):
                with m.If(self.strobe_in):
                    m.d.sync += xr_wr.data.eq(i_cooked)
                    m.d.sync += xi_wr.data.eq(q_cooked)
                    m.d.sync += xr_wr.addr.eq(revidx)
                    m.d.sync += xi_wr.addr.eq(revidx)
                    m.next = "WINDOW_WRITE"

            with m.State("WINDOW_WRITE"):
                m.d.sync += xr_wr.en.eq(1)
                m.d.sync += xi_wr.en.eq(1)
                m.d.sync += idx.eq(idx+1)
                m.next = "WINDOW"

            with m.State("FFTLOOP"):
                m.d.sync += xr_wr.en.eq(0)
                m.d.sync += xi_wr.en.eq(0)
                m.d.sync += yr_wr.en.eq(0)
                m.d.sync += yi_wr.en.eq(0)
                with m.If(idx >= pts):
                    m.d.sync += idx.eq(0)
                    m.d.sync += mask.eq(mask>>1)
                    m.d.sync += stage.eq(stage+1)

                with m.If(stage >= N):
                    m.d.sync += idx.eq(0)
                    m.next = "OUTPUT"
                with m.Else():
                    m.next = "ADDRB"

            with m.State("ADDRB"):
                with m.If(stage & 1):
                    m.d.sync += yr_rd.addr.eq(2*idx+1)
                    m.d.sync += yi_rd.addr.eq(2*idx+1)
                with m.Else():
                    m.d.sync += xr_rd.addr.eq(2*idx+1)
                    m.d.sync += xi_rd.addr.eq(2*idx+1)

                m.next = "ADDRB_LATCHED"
           
            with m.State("ADDRB_LATCHED"):
                m.next = "READB"
           
            with m.State("READB"):
                with m.If(stage & 1):
                    m.d.sync += br.eq(yr_rd.data)
                    m.d.sync += bi.eq(yi_rd.data)
                with m.Else():
                    m.d.sync += br.eq(xr_rd.data)
                    m.d.sync += bi.eq(xi_rd.data)

                m.next = "ADDRA"
           
            with m.State("ADDRA"):
                with m.If(stage & 1):
                    m.d.sync += yr_rd.addr.eq(2*idx)
                    m.d.sync += yi_rd.addr.eq(2*idx)
                with m.Else():
                    m.d.sync += xr_rd.addr.eq(2*idx)
                    m.d.sync += xi_rd.addr.eq(2*idx)

                m.next = "ADDRA_LATCHED"
           
            with m.State("ADDRA_LATCHED"):
                m.next = "READA"
           
            with m.State("READA"):
                with m.If(stage & 1):
                    m.d.sync += ar.eq(yr_rd.data)
                    m.d.sync += ai.eq(yi_rd.data)
                with m.Else():
                    m.d.sync += ar.eq(xr_rd.data)
                    m.d.sync += ai.eq(xi_rd.data)

                m.next = "BUTTERFLY"
           
            with m.State("BUTTERFLY"):
                with m.If(stage & 1):
                    m.d.sync += xr_wr.data.eq(sr)
                    m.d.sync += xi_wr.data.eq(si)
                    m.d.sync += xr_wr.addr.eq(idx)
                    m.d.sync += xi_wr.addr.eq(idx)
                with m.Else():
                    m.d.sync += yr_wr.data.eq(sr)
                    m.d.sync += yi_wr.data.eq(si)
                    m.d.sync += yr_wr.addr.eq(idx)
                    m.d.sync += yi_wr.addr.eq(idx)

                m.next = "WRITESUM"
           
            with m.State("WRITESUM"):
                with m.If(stage & 1):
                    m.d.sync += xr_wr.en.eq(1)
                    m.d.sync += xi_wr.en.eq(1)
                with m.Else():
                    m.d.sync += yr_wr.en.eq(1)
                    m.d.sync += yi_wr.en.eq(1)

                m.next = "ADDRDIFF"
           
            with m.State("ADDRDIFF"):
                with m.If(stage & 1):
                    m.d.sync += xr_wr.en.eq(0)
                    m.d.sync += xi_wr.en.eq(0)
                    m.d.sync += xr_wr.data.eq(dr)
                    m.d.sync += xi_wr.data.eq(di)
                    m.d.sync += xr_wr.addr.eq(idx+(pts>>1))
                    m.d.sync += xi_wr.addr.eq(idx+(pts>>1))
                with m.Else():
                    m.d.sync += yr_wr.en.eq(0)
                    m.d.sync += yi_wr.en.eq(0)
                    m.d.sync += yr_wr.data.eq(dr)
                    m.d.sync += yi_wr.data.eq(di)
                    m.d.sync += yr_wr.addr.eq(idx+(pts>>1))
                    m.d.sync += yi_wr.addr.eq(idx+(pts>>1))

                m.next = "WRITEDIFF"
           
            with m.State("WRITEDIFF"):
                with m.If(stage & 1):
                    m.d.sync += xr_wr.en.eq(1)
                    m.d.sync += xi_wr.en.eq(1)
                with m.Else():
                    m.d.sync += yr_wr.en.eq(1)
                    m.d.sync += yi_wr.en.eq(1)

                m.d.sync += idx.eq(idx+1)
                m.next = "FFTLOOP"

            with m.State("OUTPUT"):
                m.d.sync += self.strobe_out.eq(0)
                with m.If(idx >= pts):
                    m.next = "DONE"
                with m.Else():
                    with m.If(N & 1):
                        m.d.sync += yr_rd.addr.eq(idx)
                        m.d.sync += yi_rd.addr.eq(idx)
                    with m.Else():
                        m.d.sync += xr_rd.addr.eq(idx)
                        m.d.sync += xi_rd.addr.eq(idx)

                    m.next = "READOUT"

            with m.State("READOUT"):
                with m.If(N & 1):
                    m.d.sync += self.out_real.eq(yr_rd.data>>N)
                    m.d.sync += self.out_imag.eq(yi_rd.data>>N)
                with m.Else():
                    m.d.sync += self.out_real.eq(xr_rd.data>>N)
                    m.d.sync += self.out_imag.eq(xi_rd.data>>N)

                m.d.sync += self.strobe_out.eq(1)
                m.d.sync += idx.eq(idx+1)
                m.next = "OUTPUT"

            with m.State("DONE"):
                m.next = "IDLE"

        return m

class FixedPointFFTTest(GatewareTestCase):
    FRAGMENT_UNDER_TEST = FixedPointFFT
    FRAGMENT_ARGUMENTS = dict(bitwidth=18, pts=256)

    @sync_test_case
    def test_fft(self):
        dut = self.dut
        PTS = 256
        LPTS = log2_int(PTS)

        I =[int(cos(2*16*pi*i/PTS) * (2**16-1)) for i in range(PTS)]
        Q =[int(sin(2*16*pi*i/PTS) * (2**16-1)) for i in range(PTS)]

        for i in range(PTS):
            yield dut.in_i.eq(I[i])
            yield dut.in_q.eq(Q[i])
            yield
            yield dut.strobe_in.eq(1)
            yield
            yield dut.strobe_in.eq(0)
            yield
            yield
            yield

        # Looks that it takes ~12 cycles for read-butterfly-write
        for _ in range(PTS*LPTS*12):
            yield


if __name__ == "__main__":

    vfft = FixedPointFFT()

    ports = [
        vfft.in_i,
        vfft.in_q,
        vfft.out_real,
        vfft.out_imag,
        vfft.strobe_in,
        vfft.strobe_out,
    ]
    main(vfft, name="FixedPointFFT", ports=ports)
