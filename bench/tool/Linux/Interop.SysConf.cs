﻿// Licensed to the .NET Foundation under one or more agreements.
// The .NET Foundation licenses this file to you under the MIT license.

#nullable enable

using System.Runtime.InteropServices;

internal static partial class Interop
{
    internal static class Sys
    {
        internal enum SysConfName
        {
            _SC_CLK_TCK = 1,
            _SC_PAGESIZE = 2
        }

        [DllImport("libSystem.Native", EntryPoint = "SystemNative_SysConf", SetLastError = true)]
        internal static extern long SysConf(SysConfName name);
    }
}