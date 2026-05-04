---
title: "Third-party licenses"
sidebar_position: 6
---

# Third-Party Notices — Datagrok Platform (Core)

The Datagrok Platform Core incorporates the following third-party open source software.
Each component is provided under its own license terms, reproduced below. All
copyright notices are retained as required by their respective licenses.

This notice covers the platform core (`datlas`, `xamgle`, `ddt`, `d4`) and the
libraries vendored under `core/server/libs/` and `core/client/libs/`.

> Note: Each plugin has its own THIRD-PARTY-LICENSES.md file that is distributed along
with the plugin.

> Maintainers: regenerate this file when `pubspec.yaml` / `pubspec.lock` change in
> any of the four core packages. The companion analysis lives at
> `core/docs/deps/LICENSES.md`.

---

## Table of Contents

1. [Components under the BSD 3-Clause License](#1-bsd-3-clause)
2. [Components under the BSD 2-Clause License](#2-bsd-2-clause)
3. [Components under the MIT License](#3-mit)
4. [Components under the Apache License 2.0](#4-apache-20)
5. [Full License Texts](#5-full-license-texts)
   - [BSD 3-Clause](#bsd-3-clause-license)
   - [BSD 2-Clause](#bsd-2-clause-license)
   - [MIT](#mit-license)
   - [Apache License 2.0](#apache-license-version-20)

---

## 1. BSD 3-Clause

The following components are licensed under the BSD 3-Clause License (full text in
[Section 5](#bsd-3-clause-license)).

| Component                  | Version       | Author(s)                              |
|----------------------------|---------------|----------------------------------------|
| analyzer                   | 0.30.0+2/+4   | the Dart project authors               |
| args                       | 1.5.1         | the Dart project authors               |
| async                      | 1.13.3        | the Dart project authors               |
| benchmark_harness          | 1.0.4         | the Dart project authors               |
| browser                    | 0.10.0+3      | the Dart project authors               |
| build_runner               | 0.5.0         | the Dart project authors               |
| collection                 | 1.14.6        | the Dart project authors               |
| cookie                     | 0.0.4         | the Dart project authors               |
| crypto                     | 2.0.2+1       | the Dart project authors               |
| gcloud                     | 0.4.0+1       | the Dart project authors               |
| git                        | 0.5.1         | Kevin Moore                            |
| googleapis                 | 0.50.1        | the Dart project authors               |
| googleapis_auth            | 0.2.3+6       | the Dart project authors               |
| http                       | 0.11.3+16     | the Dart project authors               |
| intl                       | 0.15.6        | the Dart project authors               |
| js                         | 0.6.1+1       | the Dart project authors               |
| markdown                   | 1.1.1         | the Dart project authors               |
| meta                       | 1.1.8         | the Dart project authors               |
| observable                 | 0.21.0        | the Dart project authors               |
| observe                    | 0.8.10+4      | the Dart project authors               |
| path                       | 1.5.1         | the Dart project authors               |
| postgresql                 | (vendored)    | Greg Lowe                              |
| property_grid              | 1.2.8         | Ali Akbar Vathi and contributors       |
| pub_semver                 | 1.4.2         | the Dart project authors               |
| rational                   | 0.1.11        | Alexandre Ardhuin                      |
| shelf                      | 0.6.8         | the Dart project authors               |
| shelf_web_socket           | (vendored)    | the Dart project authors               |
| source_map_stack_trace     | 1.1.5         | the Dart project authors               |
| stack_trace                | 1.9.3         | the Dart project authors               |
| test                       | 0.12.24       | the Dart project authors               |
| vector_math                | 2.0.7         | Google Inc.                            |
| watcher                    | 0.9.7+6       | the Dart project authors               |
| web_socket_channel         | 1.0.6         | the Dart project authors               |

The `dock_spawn` library (vendored, version 1.0.3) by Ali Akbar Vathi is also
distributed under a BSD-style license; see `core/client/libs/dock_spawn/` for
the upstream notice.

## 2. BSD 2-Clause

The following components are licensed under the BSD 2-Clause License (full text in
[Section 5](#bsd-2-clause-license)).

| Component         | Version      | Author(s)            |
|-------------------|--------------|----------------------|
| asn1lib           | 0.4.3        | Warren Strange       |
| http_exception    | 0.1.0        | Anders Holmgren      |
| node_preamble     | 1.4.10       | Natalie Weizenbaum   |
| shelf_bind        | (vendored)   | Anders Holmgren      |
| shelf_rest        | (vendored)   | Anders Holmgren      |
| shelf_route       | (vendored)   | Anders Holmgren      |
| synchronized      | 1.5.0+1      | Alex Tekartik        |
| twitter           | 0.5.2        | Adam Singer          |

## 3. MIT

The following components are licensed under the MIT License (full text in
[Section 5](#mit-license)).

| Component         | Version       | Author(s)                                                                                |
|-------------------|---------------|------------------------------------------------------------------------------------------|
| archive           | 1.0.33        | Brendan Duncan                                                                           |
| cipher            | 0.7.1         | Iván Zaera Avellón                                                                       |
| codemirror        | 0.4.1+5.13.4  | Devon Carew (wraps CodeMirror by Marijn Haverbeke and others, MIT)                       |
| css_animation     | 0.1.8         | Devon Carew                                                                              |
| dart_amqp         | 0.3.1         | Achilleas Anagnostopoulos                                                                |
| edit_distance     | 0.3.0         | Aaron Wolen                                                                              |
| fcm_push          | 1.2.5         | Salman A. Kagzi                                                                          |
| frappe            | 0.4.0+6       | Brian Ferris                                                                             |
| image             | 1.1.33        | Brendan Duncan                                                                           |
| mailer            | 1.2.3         | the mailer contributors                                                                  |
| petitparser       | 1.7.7         | Lukas Renggli                                                                            |
| pointycastle      | 1.0.2         | The Legion of the Bouncy Castle (Australia) Pty. Ltd. and the PointyCastle project authors |
| sasl_scram        | (vendored)    | mongo-dart                                                                               |
| saslprep          | (vendored)    | Richard Burkhardt                                                                        |
| svg_pan_zoom      | 0.0.4         | Andre Tampubolon (wraps svg-pan-zoom by Andrea Leofreddi, BSD-2)                         |
| unorm_dart        | 0.3.0         | Yasuhiro Shimizu                                                                         |
| uuid              | 1.0.3         | Yulian Kuncheff                                                                          |
| xml               | 3.0.1         | Lukas Renggli                                                                            |
| xmlstream         | 0.11.1        | Stuart Sierra (Dart port)                                                                |

## 4. Apache 2.0

The following components are licensed under the Apache License, Version 2.0 (full
text in [Section 5](#apache-license-version-20)).

| Component                    | Version        | Author(s)                                            |
|------------------------------|----------------|------------------------------------------------------|
| dart_to_js_script_rewriter   | 1.0.3          | Anders Holmgren                                      |
| dartdap                      | 0.2.2          | Hoylen Sue                                           |
| google_maps                  | 3.2.2          | Alexandre Ardhuin                                    |
| junitreport                  | 0.3.3          | the junitreport authors                              |
| platform_detect              | 1.3.5          | Workiva Inc.                                         |
| quiver                       | 0.28.2 / 1.0.0 | Google Inc. and the Quiver project authors           |
| quiver_hashcode              | 1.0.0          | Google Inc.                                          |
| stream_transform             | 0.0.9          | the Dart project authors                             |

---

## 5. Full License Texts

### BSD 3-Clause License

```
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

 3. Neither the name of the copyright holder nor the names of its contributors
    may be used to endorse or promote products derived from this software
    without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
```

### BSD 2-Clause License

```
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
```

### MIT License

```
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
```

### Apache License, Version 2.0

```
                                 Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/

   TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION

   1. Definitions.

      "License" shall mean the terms and conditions for use, reproduction,
      and distribution as defined by Sections 1 through 9 of this document.

      "Licensor" shall mean the copyright owner or entity authorized by
      the copyright owner that is granting the License.

      "Legal Entity" shall mean the union of the acting entity and all
      other entities that control, are controlled by, or are under common
      control with that entity. For the purposes of this definition,
      "control" means (i) the power, direct or indirect, to cause the
      direction or management of such entity, whether by contract or
      otherwise, or (ii) ownership of fifty percent (50%) or more of the
      outstanding shares, or (iii) beneficial ownership of such entity.

      "You" (or "Your") shall mean an individual or Legal Entity
      exercising permissions granted by this License.

      "Source" form shall mean the preferred form for making modifications,
      including but not limited to software source code, documentation
      source, and configuration files.

      "Object" form shall mean any form resulting from mechanical
      transformation or translation of a Source form, including but
      not limited to compiled object code, generated documentation,
      and conversions to other media types.

      "Work" shall mean the work of authorship, whether in Source or
      Object form, made available under the License, as indicated by a
      copyright notice that is included in or attached to the work
      (an example is provided in the Appendix below).

      "Derivative Works" shall mean any work, whether in Source or Object
      form, that is based on (or derived from) the Work and for which the
      editorial revisions, annotations, elaborations, or other modifications
      represent, as a whole, an original work of authorship. For the purposes
      of this License, Derivative Works shall not include works that remain
      separable from, or merely link (or bind by name) to the interfaces of,
      the Work and Derivative Works thereof.

      "Contribution" shall mean any work of authorship, including
      the original version of the Work and any modifications or additions
      to that Work or Derivative Works thereof, that is intentionally
      submitted to Licensor for inclusion in the Work by the copyright owner
      or by an individual or Legal Entity authorized to submit on behalf of
      the copyright owner. For the purposes of this definition, "submitted"
      means any form of electronic, verbal, or written communication sent
      to the Licensor or its representatives, including but not limited to
      communication on electronic mailing lists, source code control systems,
      and issue tracking systems that are managed by, or on behalf of, the
      Licensor for the purpose of discussing and improving the Work, but
      excluding communication that is conspicuously marked or otherwise
      designated in writing by the copyright owner as "Not a Contribution."

      "Contributor" shall mean Licensor and any individual or Legal Entity
      on behalf of whom a Contribution has been received by Licensor and
      subsequently incorporated within the Work.

   2. Grant of Copyright License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      copyright license to reproduce, prepare Derivative Works of,
      publicly display, publicly perform, sublicense, and distribute the
      Work and such Derivative Works in Source or Object form.

   3. Grant of Patent License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      (except as stated in this section) patent license to make, have made,
      use, offer to sell, sell, import, and otherwise transfer the Work,
      where such license applies only to those patent claims licensable
      by such Contributor that are necessarily infringed by their
      Contribution(s) alone or by combination of their Contribution(s)
      with the Work to which such Contribution(s) was submitted. If You
      institute patent litigation against any entity (including a
      cross-claim or counterclaim in a lawsuit) alleging that the Work
      or a Contribution incorporated within the Work constitutes direct
      or contributory patent infringement, then any patent licenses
      granted to You under this License for that Work shall terminate
      as of the date such litigation is filed.

   4. Redistribution. You may reproduce and distribute copies of the
      Work or Derivative Works thereof in any medium, with or without
      modifications, and in Source or Object form, provided that You
      meet the following conditions:

      (a) You must give any other recipients of the Work or
          Derivative Works a copy of this License; and

      (b) You must cause any modified files to carry prominent notices
          stating that You changed the files; and

      (c) You must retain, in the Source form of any Derivative Works
          that You distribute, all copyright, patent, trademark, and
          attribution notices from the Source form of the Work,
          excluding those notices that do not pertain to any part of
          the Derivative Works; and

      (d) If the Work includes a "NOTICE" text file as part of its
          distribution, then any Derivative Works that You distribute must
          include a readable copy of the attribution notices contained
          within such NOTICE file, excluding those notices that do not
          pertain to any part of the Derivative Works, in at least one
          of the following places: within a NOTICE text file distributed
          as part of the Derivative Works; within the Source form or
          documentation, if provided along with the Derivative Works; or,
          within a display generated by the Derivative Works, if and
          wherever such third-party notices normally appear. The contents
          of the NOTICE file are for informational purposes only and
          do not modify the License. You may add Your own attribution
          notices within Derivative Works that You distribute, alongside
          or as an addendum to the NOTICE text from the Work, provided
          that such additional attribution notices cannot be construed
          as modifying the License.

      You may add Your own copyright statement to Your modifications and
      may provide additional or different license terms and conditions
      for use, reproduction, or distribution of Your modifications, or
      for any such Derivative Works as a whole, provided Your use,
      reproduction, and distribution of the Work otherwise complies with
      the conditions stated in this License.

   5. Submission of Contributions. Unless You explicitly state otherwise,
      any Contribution intentionally submitted for inclusion in the Work
      by You to the Licensor shall be under the terms and conditions of
      this License, without any additional terms or conditions.
      Notwithstanding the above, nothing herein shall supersede or modify
      the terms of any separate license agreement you may have executed
      with Licensor regarding such Contributions.

   6. Trademarks. This License does not grant permission to use the trade
      names, trademarks, service marks, or product names of the Licensor,
      except as required for describing the origin of the Work and
      reproducing the content of the NOTICE file.

   7. Disclaimer of Warranty. Unless required by applicable law or
      agreed to in writing, Licensor provides the Work (and each
      Contributor provides its Contributions) on an "AS IS" BASIS,
      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
      implied, including, without limitation, any warranties or conditions
      of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A
      PARTICULAR PURPOSE. You are solely responsible for determining the
      appropriateness of using or redistributing the Work and assume any
      risks associated with Your exercise of permissions under this License.

   8. Limitation of Liability. In no event and under no legal theory,
      whether in tort (including negligence), contract, or otherwise,
      unless required by applicable law (such as deliberate and grossly
      negligent acts) or agreed to in writing, shall any Contributor be
      liable to You for damages, including any direct, indirect, special,
      incidental, or consequential damages of any character arising as a
      result of this License or out of the use or inability to use the
      Work (including but not limited to damages for loss of goodwill,
      work stoppage, computer failure or malfunction, or any and all
      other commercial damages or losses), even if such Contributor
      has been advised of the possibility of such damages.

   9. Accepting Warranty or Support. While redistributing the Work or
      Derivative Works thereof, You may choose to offer, and charge a
      fee for, acceptance of support, warranty, indemnity, or other
      liability obligations and/or rights consistent with this License.
      However, in accepting such obligations, You may act only on Your
      own behalf and on Your sole responsibility, not on behalf of any
      other Contributor, and only if You agree to indemnify, defend,
      and hold each Contributor harmless for any liability incurred by,
      or claims asserted against, such Contributor by reason of your
      accepting any such warranty or support.

   END OF TERMS AND CONDITIONS
```

---

*Last updated: 2026-05-04. If you believe a notice or attribution is missing or
incorrect, please contact info@datagrok.ai.*
