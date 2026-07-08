# Profile Settings

Users can change the information in their profiles. Such as profile photo, password, first and last names.

## Testing scenarios

1. Edit name. Profile information is changed by the entered data and persists.

> Automated by `profile-settings-spec.ts` (name edit via `dapi.users.save`,
> verified across a re-fetch, original name restored afterwards).
>
> Profile-photo upload (multiple formats + invalid-format rejection) and the
> "Change password..." negative validation (mismatched / wrong current
> password) are manual: see `profile-settings-ui.md`.

---
{
  "order": 7
}